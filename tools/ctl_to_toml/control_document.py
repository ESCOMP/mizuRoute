"""CIME-independent parsing for flat mizuRoute settings documents."""

from dataclasses import dataclass
from datetime import date, datetime, time
import json
import re
import tomllib


_NAME = re.compile(r"[A-Za-z0-9_-]+\Z")
_INTEGER = re.compile(r"[+-]?[0-9]+\Z")
_FLOAT = re.compile(
    r"[+-]?(?:(?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?"
    r"|[0-9]+[eEdD][+-]?[0-9]+)\Z"
)


def parse_legacy_value(value):
    """Infer a scalar type; explicitly quoted legacy values remain strings."""
    if value.startswith('"'):
        parsed = tomllib.loads("value = " + value)["value"]
        if not isinstance(parsed, str):
            raise ValueError("Expected a quoted string")
        return parsed
    logical = value.lower()
    if logical in ("t", "true", ".true."):
        return True
    if logical in ("f", "false", ".false."):
        return False
    if logical in ("inf", "+inf", "-inf", "nan", "+nan", "-nan"):
        return float(logical)
    if _INTEGER.fullmatch(value):
        return int(value)
    if _FLOAT.fullmatch(value):
        return float(value.lower().replace("d", "e"))
    return value


def _split_comment(text, marker, toml=False):
    """Find a comment delimiter outside quoted strings (including TOML triples)."""
    quote = None
    index = 0
    while index < len(text):
        if quote:
            if quote.startswith('"') and text[index] == "\\":
                index += 2
                continue
            if text.startswith(quote, index):
                index += len(quote)
                quote = None
                continue
        else:
            if text[index] == marker:
                return text[:index], text[index + 1:]
            if (text[index] == '"' and (toml or not text[:index].strip())) or (toml and text[index] == "'"):
                char = text[index]
                quote = char * 3 if toml and text.startswith(char * 3, index) else char
                index += len(quote)
                continue
        index += 1
    return text, None


@dataclass
class Entry:
    """A blank line, comment, or scalar setting in source order."""

    key: str = None
    value: object = None
    comment: str = None
    indent: str = ""


class ControlDocument:
    """A settings document retaining comment text and its position among keys."""

    def __init__(self, entries):
        self.entries = entries

    @classmethod
    def parse(cls, text, suffix):
        if suffix == ".control":
            entries = cls._parse_control(text)
        elif suffix == ".toml":
            entries = cls._parse_toml(text)
        else:
            raise ValueError("Expected a .control or .toml extension")
        if not any(entry.key is not None for entry in entries):
            raise ValueError("No settings were found in the input file")
        return cls(entries)

    @staticmethod
    def _parse_control(text):
        entries = []
        seen = set()
        for number, line in enumerate(text.splitlines(), 1):
            indent = line[:len(line) - len(line.lstrip())]
            if not line.strip():
                entries.append(Entry(indent=line))
                continue
            if line.lstrip().startswith("!"):
                entries.append(Entry(comment=line.lstrip()[1:], indent=indent))
                continue
            match = re.fullmatch(r"\s*<([^>]+)>\s+(.*)", line)
            if not match:
                raise ValueError(f"Line {number}: expected <setting> value ! comment")
            key = match[1].strip()
            if not _NAME.fullmatch(key) or key in seen:
                raise ValueError(f"Line {number}: invalid or duplicate setting {key!r}")
            raw, comment = _split_comment(match[2], "!")
            raw = raw.strip()
            # Move leading zeros to the end so integer conversion keeps
            # routing method 0 (runoff accumulation) enabled.
            if key == "route_opt" and re.fullmatch(r"\+?0[0-9]+", raw):
                digits = raw.lstrip("+")
                nonzero_prefix = digits.lstrip("0")
                raw = nonzero_prefix + "0" * (len(digits) - len(nonzero_prefix))
            try:
                value = parse_legacy_value(raw)
            except (ValueError, KeyError) as error:
                raise ValueError(f"Line {number}: {error}") from error
            entries.append(Entry(key, value, comment, indent))
            seen.add(key)
        return entries

    @staticmethod
    def _parse_toml(text):
        # Validate the entire document first, including duplicate keys.
        parsed = tomllib.loads(text)
        entries = []
        lines = iter(enumerate(text.splitlines(), 1))
        for number, line in lines:
            indent = line[:len(line) - len(line.lstrip())]
            if not line.strip():
                entries.append(Entry(indent=line))
                continue
            if line.lstrip().startswith("#"):
                entries.append(Entry(comment=line.lstrip()[1:], indent=indent))
                continue
            statement = line
            while True:
                try:
                    setting = tomllib.loads(statement)
                    break
                except tomllib.TOMLDecodeError:
                    try:
                        _, continuation = next(lines)
                    except StopIteration as error:
                        raise ValueError(f"Line {number}: incomplete TOML setting") from error
                    statement += "\n" + continuation
            if len(setting) != 1:
                raise ValueError(f"Line {number}: expected a flat scalar setting")
            key, value = next(iter(setting.items()))
            if not _NAME.fullmatch(key) or isinstance(value, (dict, list)):
                raise ValueError(f"Line {number}: only flat scalar settings are supported")
            if isinstance(value, (datetime, date, time)):
                value = value.isoformat(sep=" ") if isinstance(value, datetime) else value.isoformat()
            _, comment = _split_comment(statement, "#", toml=True)
            entries.append(Entry(key, value, comment, indent))
        if len(parsed) != sum(entry.key is not None for entry in entries):
            raise ValueError("Only flat scalar settings are supported")
        return entries

    def render(self, suffix):
        """Serialize values and translate comment delimiters without losing text."""
        marker = "#" if suffix == ".toml" else "!"
        lines = []
        for entry in self.entries:
            if entry.key is None:
                lines.append(entry.indent + (marker + entry.comment if entry.comment is not None else ""))
                continue
            value = entry.value
            if isinstance(value, bool):
                formatted = ("true" if value else "false") if suffix == ".toml" else ("T" if value else "F")
            elif isinstance(value, (int, float)):
                formatted = str(value)
            elif suffix == ".toml":
                formatted = json.dumps(value, ensure_ascii=True)
            else:
                formatted = value
                # Legacy output cannot escape these characters without quotes.
                # Refuse values that would be truncated or reinterpreted on read.
                if any(char in formatted for char in "!\r\n") or formatted.startswith('"') or formatted != formatted.strip():
                    raise ValueError(f"Setting {entry.key!r} cannot be represented as an unquoted legacy string")
            if suffix == ".toml":
                line = f"{entry.indent}{entry.key} = {formatted}"
            else:
                line = f"{entry.indent}<{entry.key}> {formatted}"
            if entry.comment is not None:
                line += " " + marker + entry.comment
            elif suffix == ".control":
                # The historical reader requires a trailing comment delimiter.
                line += " !"
            lines.append(line)
        output = "\n".join(lines) + "\n"
        if suffix == ".toml":
            tomllib.loads(output)
        return output
