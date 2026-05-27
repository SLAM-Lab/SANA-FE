"""Generate models.rst from sanafe introspection.

Run after any change to a model's supported_attributes map:
    python docs/make_model_descriptions.py

The generated file should be committed alongside the model changes that
produced it.
"""
import sanafe
import os

OUTPUT_PATH = "models.rst"

PAGE_TEMPLATE = """\
Models
======

This page documents every pipeline model available in SANA-FE. It begins with
the framework attributes shared by all models, followed by a section for each
model listing its model-specific attributes.


Framework attributes
--------------------

These attributes are accepted by every pipeline model in SANA-FE. They control
shared behavior such as energy/latency reporting, hardware unit selection, and
update scheduling, and apply equally to built-in and plugin models.

{framework_table}

{model_sections}
"""

MODEL_SECTION_TEMPLATE = """\
{title}
{underline}

{model_table}
"""

TABLE_TEMPLATE = """\
.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
{rows}
"""


def format_table(attrs: dict) -> str:
    if not attrs:
        return "*(none)*"
    rows = "\n".join(
        f"   * - ``{name}``\n     - {desc or '*(no description)*'}"
        for name, desc in sorted(attrs.items())
    )
    return TABLE_TEMPLATE.format(rows=rows)


def generate_model_section(model_name: str) -> str:
    info = sanafe.model_attributes[model_name]
    return MODEL_SECTION_TEMPLATE.format(
        title=model_name,
        underline="-" * len(model_name),
        model_table=format_table(info),
    )


def generate_page() -> str:
    model_sections = "\n\n".join(
        generate_model_section(name) for name in sanafe.model_attributes
    )
    return PAGE_TEMPLATE.format(
        framework_table=format_table(sanafe.framework_attributes),
        model_sections=model_sections,
    )


def main():
    out_dir = os.path.dirname(OUTPUT_PATH)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(OUTPUT_PATH, "w") as out_file:
        out_file.write(generate_page())
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
