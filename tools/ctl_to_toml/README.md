# mizuRoute settings converter

Convert mizuRoute settings between legacy `.control` and `.toml` files. Requires Python 3.11 or later, or Python 3.9–3.10 with `tomli` installed (`python3 -m pip install tomli`). This tool does not require CIME python libraries.

From the mizuRoute sandbox root:

```bash
tools/ctl_to_toml/settings_converter path/to/settings.control
tools/ctl_to_toml/settings_converter --verbose path/to/settings.toml
tools/ctl_to_toml/settings_converter --help
```

You can also run the script with `python3`, or add `tools/ctl_to_toml` to your `PATH` to use `settings_converter` directly. It works from any working directory within the sandbox checkout or outside it.

The output is saved beside the input with the opposite extension (`.control` -> `.toml` or `.toml` -> `.control`). Existing files are never overwritten: the converter tries `settings.toml`, `settings_YYYYMMDD.toml`, `settings_YYYYMMDD_1.toml`, and increasing numbers. Dates use the machine's local date. The same rule applies to `.control` output. Successful conversion is silent by default; `-v`/`--verbose` logs progress and the output path to stderr. Errors always go to stderr with a nonzero exit status.

Standalone comments, section headings, blank lines, inline comment text, and setting order are preserved. Comment delimiters change between `!` and `#`; alignment and numeric spelling may change. Integers remain integers and floating-point values remain floats, including scientific notation. Legacy `T`/`F` become TOML `true`/`false`; legacy logical spellings `true`, `false`, `.true.`, and `.false.` are also accepted. Date and time values become strings.

Double-quoted legacy values remain strings when converted to TOML. For example, the quotes in this legacy setting keep `0012` as text rather than the number `12`:

```text
<case_name> "0012" ! simulation name
```

The TOML output is:

```toml
case_name = "0012" # simulation name
```

When writing TOML strings, the converter adds surrounding quotes and uses escape sequences for special characters, such as backslashes and quotation marks, to preserve the text. For example, an unquoted legacy path:

```text
<input_dir> C:\data ! input directory
```

becomes:

```toml
input_dir = "C:\\data" # input directory
```

In the TOML file, `\\` represents one actual backslash, so the path still means `C:\data`.

TOML strings lose their surrounding quotes in legacy output. Numeric-looking strings may therefore become numbers on a subsequent conversion; strings such as `"T"` may similarly become logical values. Strings containing `!`, line breaks, leading/trailing whitespace, or an initial double quote cannot be represented safely without quotes in legacy format and produce an error. Nested TOML tables and arrays are unsupported because mizuRoute settings are flat scalars. Invalid, duplicate, or empty settings files are rejected before output is created.

The Command Line Interface (CLI) uses the local parser in `tools/ctl_to_toml/control_document.py`. The converter directory can be copied outside the repository and run independently; keep `settings_converter` and `control_document.py` together. No additional TOML writer or comment-parser dependency is required. Repository samples are used only as read-only test fixtures.

Run the regression tests from the sandbox root:

```bash
python3 -m unittest discover -s tools/ctl_to_toml/testing -p 'test_*.py' -v
```
