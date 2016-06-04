#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import unicode_literals

import re

from mathics.core.parser.errors import ScanError, IncompleteSyntaxError
from mathics.core.parser.prescanner import prescan
from mathics.core.characters import letters, letterlikes, named_characters


# special patterns
number_pattern = r'''
( (?# Two possible forms depending on whether base is specified)
    (\d+\^\^([a-zA-Z0-9]+\.?[a-zA-Z0-9]*|[a-zA-Z0-9]*\.?[a-zA-Z0-9]+))
    | (\d+\.?\d*|\d*\.?\d+)
)
(``?(\+|-)?(\d+\.?\d*|\d*\.?\d+)|`)?        (?# Precision / Accuracy)
(\*\^(\+|-)?\d+)?                           (?# Exponent)
'''
base_symbol_pattern = r'((?![0-9])([0-9${0}{1}])+)'.format(letters, letterlikes)
full_symbol_pattern = r'(`?{0}(`{0})*)'.format(base_symbol_pattern)
pattern_pattern = r'{0}?_(\.|(__?)?{0}?)?'.format(full_symbol_pattern)

tokens = [
    ('Number', number_pattern),
    ('String', r'"'),
    ('Pattern', pattern_pattern),
    ('Symbol', full_symbol_pattern),
    ('SlotSequence', r'\#\#\d*'),
    ('Slot', r'\#\d*'),
    ('Out', r'\%(\%+|\d+)?'),
    ('PutAppend', r'\>\>\>'),
    ('Put', r'\>\>'),
    ('Get', r'\<\<'),
    ('RawLeftBracket', r' \[ '),
    ('RawRightBracket', r' \] '),
    ('RawLeftBrace', r' \{ '),
    ('RawRightBrace', r' \} '),
    ('RawLeftParenthesis', r' \( '),
    ('RawRightParenthesis', r' \) '),

    ('RawComma', r' \, '),

    ('Span', r' \;\; '),

    ('MessageName', r' \:\: '),

    # # Box Constructors
    # ('InterpretedBox', r' \\\! '),
    # ('boxes_Superscript', r' \\\^ '),
    # ('boxes_Subscript', r' \\\_ '),
    # ('boxes_Overscript', r' \\\& '),
    # ('boxes_Underscript', r' \\\+ '),
    # ('boxes_Otherscript', r' \\\% '),
    # ('boxes_Fraction', r' \\\/ '),
    # ('boxes_Sqrt', r' \\\@ '),
    # ('boxes_FormBox', r' \\\` '),

    ('PatternTest', r' \? '),
    ('Increment', r' \+\+ '),
    ('Decrement', r' \-\- '),

    ('MapAll', r' \/\/\@ '),
    ('Map', r' \/\@ '),
    ('ApplyList', r' \@\@\@ '),
    ('Apply', r' \@\@ '),
    ('Prefix', r' \@ '),

    ('StringExpression', r' \~\~ '),
    ('Infix', r' \~ '),

    ('Factorial2', r' \!\! '),
    ('Factorial', r' \! '),

    ('Derivative', r' \'+ '),
    ('StringJoin', r' \<\> '),

    ('Power', r' \^ '),

    ('NonCommutativeMultiply', r' \*\* '),

    ('AddTo', r' \+\= '),
    ('SubtractFrom', r' \-\=  '),
    ('TimesBy', r' \*\= '),
    ('DivideBy', r' \/\=  '),

    ('RawDot', r' \. '),

    ('Plus', r' \+ '),
    ('Minus', r' \- '),
    ('RawBackslash', r' \\ '),

    ('Times', r' \*|\u00d7 '),

    ('Divide', r' \/|\u00f7 '),

    ('SameQ', r' \=\=\= '),
    ('UnsameQ', r' \=\!\= '),

    ('op_Equal', r' \=\= '),
    ('op_Unequal', r' \!\= '),
    ('op_GreaterEqual', r' \>\= '),
    ('op_LessEqual', r' \<\= '),
    ('Greater', r' \> '),
    ('Less', r' \< '),

    ('op_And', r' \&\& '),
    ('op_Or', r' \|\|  '),

    ('RepeatedNull', r' \.\.\. '),
    ('Repeated', r' \.\. '),
    ('Alternatives', r' \| '),

    ('RawColon', r' \: '),
    ('Condition', r' \/\; '),

    ('op_Rule', r' \-\> '),
    ('op_RuleDelayed', r' \:\> '),
    ('ReplaceRepeated', r' \/\/\. '),
    ('ReplaceAll', r' \/\. '),

    ('RawAmpersand', r' \& '),
    ('Postfix', r' \/\/ '),

    ('UpSetDelayed', r' \^\:\= '),
    ('SetDelayed', r' \:\= '),
    ('UpSet', r' \^\= '),
    ('TagSet', r' \/\: '),
    ('Unset', r' \=\. '),
    ('Set', r' \= '),

    ('Semicolon', r' \; '),

    # ('DiscreteShift', r' \uf4a3 '),
    # ('DiscreteRatio', r' \uf4a4 '),
    # ('DifferenceDelta', r' \u2206 '),
    # ('PartialD', r' \u2202 '),

    ('Cross', r' \uf4a0 '),
    ('Colon', r' \u2236 '),
    ('Transpose', r' \uf3c7 '),
    ('Conjugate', r' \uf3c8 '),
    ('ConjugateTranspose', r' \uf3c9 '),
    ('HermitianConjugate', r' \uf3ce '),
    ('Integral', r' \u222b '),
    ('DifferentialD', r' \uf74c '),
    ('Del', r' \u2207 '),
    ('Square', r' \uf520 '),
    ('SmallCircle', r' \u2218 '),
    ('CircleDot', r' \u2299 '),

    # ('Sum', r' \u2211 '),
    # ('Product', r' \u220f '),
    ('PlusMinus', r' \u00b1 '),
    ('MinusPlus', r' \u2213 '),
    ('Or', r' \u2228 '),
    ('Nor', r' \u22BD '),
    ('And', r' \u2227 '),
    ('Nand', r' \u22BC '),
    ('Xor', r' \u22BB '),
    ('Xnor', r' \uF4A2 '),
    ('Diamond', r' \u22c4 '),
    ('Wedge', r' \u22c0 '),
    ('Vee', r' \u22c1 '),
    ('CircleTimes', r' \u2297 '),
    ('CenterDot', r' \u00b7 '),
    ('Star', r' \u22c6'),
    ('VerticalTilde', r' \u2240 '),
    ('Coproduct', r' \u2210 '),
    ('Cap', r' \u2322 '),
    ('Cup', r' \u2323 '),
    ('CirclePlus', r' \u2295 '),
    ('CircleMinus', r' \u2296 '),
    ('Intersection', r' \u22c2 '),
    ('Union', r' \u22c3 '),
    ('Equal', r' \uf431 '),
    ('LongEqual', r' \uf7d9 '),
    ('NotEqual', r' \u2260 '),
    ('LessEqual', r' \u2264 '),
    ('LessSlantEqual', r' \u2a7d '),
    ('GreaterEqual', r' \u2265 '),
    ('GreaterSlantEqual', r' \u2a7e '),
    ('VerticalBar', r' \u2223 '),
    ('NotVerticalBar', r' \u2224 '),
    ('DoubleVerticalBar', r' \u2225 '),
    ('NotDoubleVerticalBar', r' \u2226 '),
    ('Element', r' \u2208 '),
    ('NotElement', r' \u2209 '),
    ('Subset', r' \u2282 '),
    ('Superset', r' \u2283 '),
    ('ForAll', r' \u2200 '),
    ('Exists', r' \u2203 '),
    ('NotExists', r' \u2204 '),
    ('Not', r' \u00AC '),
    ('Equivalent', r' \u29E6 '),
    ('Implies', r' \uF523 '),
    ('RightTee', r' \u22A2 '),
    ('DoubleRightTee', r' \u22A8 '),
    ('LeftTee', r' \u22A3 '),
    ('DoubleLeftTee', r' \u2AE4 '),
    ('SuchThat', r' \u220D '),
    ('Rule', r' \uF522 '),
    ('RuleDelayed', r' \uF51F '),
    ('VerticalSeparator', r' \uF432 '),
    ('Therefore', r' \u2234 '),
    ('Because', r' \u2235 '),
    ('Function', r' \uF4A1 '),
]


# compile tokens
tokens = [(tag, re.compile(pattern, re.VERBOSE)) for tag, pattern in tokens]

filename_pattern = re.compile(
    r'''
    (?P<quote>\"?)                              (?# Opening quotation mark)
        [a-zA-Z0-9\`/\.\\\!\-\:\_\$\*\~\?]+     (?# Literal characters)
    (?P=quote)                                  (?# Closing quotation mark)
    ''', re.VERBOSE)


class Token(object):
    def __init__(self, tag, text, pos):
        self.tag = tag
        self.text = prescan(text)
        self.pos = pos


class Tokeniser(object):
    def __init__(self, code):
        self.pos = 0
        self.code = code

    def next(self, tokens=tokens):
        'return next token'
        self.skip_blank()
        if self.pos >= len(self.code):
            return Token('END', '', len(self.code) - 1)
        for tag, pattern in tokens:
            match = pattern.match(self.code, self.pos)
            if match is not None:
                # look for custom tokenisation rule
                override = getattr(self, 't_' + tag, None)
                if override is not None:
                    return override(match)
                else:
                    text = match.group(0)
                    self.pos = match.end(0)
                    return Token(tag, text, self.pos)
        raise ScanError(self.pos)

    def next_filename(self):
        'return next filename token'
        return self.next(tokens=[('filename', filename_pattern)])

    def skip_blank(self):
        'skip whitespace and comments'
        comment = []   # start positions of comments
        while self.pos < len(self.code):
            if comment:
                if self.code.startswith('(*', self.pos):
                    comment.append(self.pos)
                    self.pos += 2
                elif self.code.startswith('*)', self.pos):
                    comment.pop()
                    self.pos += 2
                else:
                    self.pos += 1
            elif self.code.startswith('(*', self.pos):
                comment.append(self.pos)
                self.pos += 2
            elif self.code[self.pos] in ' \r\n\t':
                self.pos += 1
            else:
                break
        if comment:
            self.pos = comment[0]
            raise IncompleteSyntaxError()

    def t_String(self, match):
        start, end = self.pos, None
        self.pos += 1   # skip opening '"'
        while self.pos < len(self.code):
            c = self.code[self.pos]
            if c == '"':
                self.pos += 1
                end = self.pos
                break
            elif c == '\\':
                self.pos += 2
            else:
                self.pos += 1
        if end is None:
            # reached end while still inside string
            self.pos = start - 1
            raise IncompleteSyntaxError(self.pos)
        return Token('String', self.code[start:end], start)
