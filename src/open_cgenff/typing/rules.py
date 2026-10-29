"""Parser for the typing rule language (P-I; SPEC 5.5).

    cat <NAME>                      category; ends at `end` or at the next `cat`
    typ <TYPE> : <items>            assign TYPE when all conditions hold
    sub <CAT>  : <items>            continue in category CAT
    def <NAME> ( <items> )          macro (extension): `@NAME` anywhere expands to <items>
    # comment                       to end of line; a rule continues on the next line while parentheses are open

Conditions (P-I Table 1): el X, elha, elos, nb n, rings n, ring3/ring2/arom/ring23/ring [s], self, ne (g) (g)...,
bo n, inring, ! (g), or (g) (g)... Optional actions (top level of a rule only): charge n, impr, warn "msg",
err "msg", altnum.

Extensions (not in P-I; SPEC 5.5 lists them):
  nh n          exactly n hydrogen neighbours            deg n        exactly n neighbours
  fc n          resonance formal charge equals n (SPEC 5.3)
  inring [class] [s]   the last bond traversed by `ne` lies in a ring of that class and/or size
  type PATTERN  the atom's own type (typed on demand) matches a glob pattern, e.g. `type CG2R6*`
  rmember class [s] (g)   some atom of one of the atom's rings of that class/size satisfies g
  fusion        bridgehead whose rings (SPEC 5.2 ring list) share a bond or bridge, i.e. not a spiro atom
Ring keywords without a size match any size. Ring class keywords: ring3 (all sp3), ring2 (all sp2), arom, ring23
(mixed), ring (any).
"""
import re
from dataclasses import dataclass, field

RING_CLASSES = {"ring3": "sp3", "ring2": "sp2", "arom": "arom", "ring23": "mixed", "ring": None}
_INT_ARG = {"nb", "rings", "bo", "nh", "deg", "fc"}
_TOKEN = re.compile(r'"[^"]*"|[()]|[^\s()"]+')


class RuleSyntaxError(ValueError):
    pass


@dataclass
class Rule:
    action: str                       # 'typ' | 'sub'
    target: str                       # type or category name
    conditions: list                  # condition tuples, see _parse_items
    charge: int = None
    impr: bool = False
    warn: list = field(default_factory=list)
    err: list = field(default_factory=list)
    altnum: bool = False
    line: int = 0
    category: str = ""

    def __str__(self):
        return f"{self.category}:{self.line} {self.action} {self.target}"


@dataclass
class RuleSet:
    categories: dict                  # name -> [Rule], in file order
    source: str = ""

    @classmethod
    def from_file(cls, path):
        with open(path) as fh:
            return parse(fh.read(), str(path))

    def types(self):
        return sorted({r.target for rules in self.categories.values() for r in rules if r.action == "typ"})


def _logical_lines(text):
    """(line number, tokens) per rule, joining lines while parentheses are open."""
    buf, depth, start = [], 0, 0
    for n, raw in enumerate(text.splitlines(), 1):
        toks = []
        for t in _TOKEN.findall(raw):                   # quoted messages may contain '#'
            if t.startswith("#"):
                break
            if "#" in t and not t.startswith('"'):
                toks.append(t.split("#", 1)[0])
                break
            toks.append(t)
        if not toks:
            continue
        if not buf:
            start = n
        buf += toks
        depth += toks.count("(") - toks.count(")")
        if depth < 0:
            raise RuleSyntaxError(f"line {n}: unbalanced ')'")
        if depth == 0:
            yield start, buf
            buf = []
    if buf:
        raise RuleSyntaxError(f"line {start}: unclosed '('")


class _Parser:
    def __init__(self, toks, line, macros):
        self.toks, self.k, self.line, self.macros = toks, 0, line, macros

    def error(self, msg):
        raise RuleSyntaxError(f"line {self.line}: {msg}")

    def peek(self):
        return self.toks[self.k] if self.k < len(self.toks) else None

    def take(self):
        t = self.peek()
        if t is None:
            self.error("unexpected end of rule")
        self.k += 1
        return t

    def take_int(self, kw):
        t = self.take()
        try:
            return int(t)
        except ValueError:
            self.error(f"'{kw}' needs an integer, got {t!r}")

    def group(self):
        if self.take() != "(":
            self.error("expected '('")
        items, opts = self.items(until=")")
        if opts:
            self.error("actions (charge, impr, warn, err, altnum) are only allowed at the top level of a rule")
        self.take()
        return items

    def groups(self):
        gs = []
        while self.peek() == "(":
            gs.append(self.group())
        if not gs:
            self.error("expected at least one '( ... )' group")
        return gs

    def optional_size(self):
        t = self.peek()
        if t is not None and t.lstrip("-").isdigit():
            return int(self.take())
        return None

    def items(self, until=None):
        conds, opts = [], {}
        while self.peek() is not None and self.peek() != until:
            kw = self.take()
            if kw == "el":
                conds.append(("el", self.take()))
            elif kw in ("elha", "elos", "self", "fusion"):
                conds.append((kw,))
            elif kw in _INT_ARG:
                conds.append((kw, self.take_int(kw)))
            elif kw in RING_CLASSES:
                conds.append(("ringkw", RING_CLASSES[kw], self.optional_size()))
            elif kw == "inring":
                cls = RING_CLASSES[self.take()] if self.peek() in RING_CLASSES else None
                conds.append(("inring", cls, self.optional_size()))
            elif kw == "rmember":
                c = self.take()
                if c not in RING_CLASSES:
                    self.error(f"'rmember' needs a ring class, got {c!r}")
                size = self.optional_size()
                conds.append(("rmember", RING_CLASSES[c], size, self.group()))
            elif kw == "type":
                conds.append(("type", self.take()))
            elif kw == "ne":
                conds.append(("ne", self.groups()))
            elif kw == "or":
                conds.append(("or", self.groups()))
            elif kw == "!":
                conds.append(("not", self.group()))
            elif kw.startswith("@"):
                if kw[1:] not in self.macros:
                    self.error(f"undefined macro {kw}")
                conds.extend(self.macros[kw[1:]])
            elif kw == "charge":
                opts["charge"] = self.take_int(kw)
            elif kw in ("impr", "altnum"):
                opts[kw] = True
            elif kw in ("warn", "err"):
                msg = self.take()
                if not (msg.startswith('"') and msg.endswith('"')):
                    self.error(f"'{kw}' needs a quoted message")
                opts.setdefault(kw, []).append(msg[1:-1])
            else:
                self.error(f"unknown keyword {kw!r}")
        return conds, opts


def parse(text, source=""):
    categories, macros, cat = {}, {}, None
    pending_subs = []
    for line, toks in _logical_lines(text):
        head = toks[0]
        if head == "cat":
            if len(toks) != 2:
                raise RuleSyntaxError(f"line {line}: 'cat <NAME>'")
            cat = toks[1]
            if cat in categories:
                raise RuleSyntaxError(f"line {line}: category {cat} defined twice")
            categories[cat] = []
        elif head == "end":
            cat = None
        elif head == "def":
            if len(toks) < 4 or toks[2] != "(" or toks[-1] != ")":
                raise RuleSyntaxError(f"line {line}: 'def <NAME> ( <conditions> )'")
            p = _Parser(toks[2:], line, macros)
            macros[toks[1]] = p.group()
            if p.peek() is not None:
                raise RuleSyntaxError(f"line {line}: text after macro body")
        elif head in ("typ", "sub"):
            if cat is None:
                raise RuleSyntaxError(f"line {line}: rule outside a category")
            if len(toks) < 3 or toks[2] != ":":
                raise RuleSyntaxError(f"line {line}: '{head} <NAME> : <conditions>'")
            conds, opts = _Parser(toks[3:], line, macros).items()
            rule = Rule(head, toks[1], conds, opts.get("charge"), opts.get("impr", False), opts.get("warn", []),
                        opts.get("err", []), opts.get("altnum", False), line, cat)
            categories[cat].append(rule)
            if head == "sub":
                pending_subs.append(rule)
        else:
            raise RuleSyntaxError(f"line {line}: expected cat/typ/sub/def/end, got {head!r}")
    if "main" not in categories:
        raise RuleSyntaxError("no 'cat main'")
    for rule in pending_subs:
        if rule.target not in categories:
            raise RuleSyntaxError(f"line {rule.line}: sub to undefined category {rule.target}")
    return RuleSet(categories, source)
