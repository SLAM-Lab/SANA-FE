"""Generate model documentation .rst files from sanafe introspection.

Run after any change to a model's supported_attributes map:
    python docs/make_model_descriptions.py

Generated files should be committed alongside the model changes that produced them.
"""
import sanafe
import os

OUTPUT_DIR = "models"

MODEL_PAGE_TEMPLATE = """\
{title}
{underline}


Model-specific attributes
-------------------------

{model_table}

Framework attributes
--------------------

These attributes are shared by all pipeline models. See
:doc:`framework_attributes` for the full reference.

{framework_table}
"""

TABLE_TEMPLATE = """\
.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
{rows}
"""

FRAMEWORK_PAGE_TEMPLATE = """\
Framework attributes
====================

These attributes are accepted by every pipeline model in SANA-FE. They control
shared behavior such as energy/latency reporting, hardware unit selection, and
update scheduling, and apply equally to built-in and plugin models.

For attributes specific to a particular model, see that model's documentation
page.

{framework_table}
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


def generate_model_page(model_name: str) -> str:
    info = sanafe.model_attributes[model_name]
    return MODEL_PAGE_TEMPLATE.format(
        title=model_name,
        underline="=" * len(model_name),
        model_table=format_table(info),
        framework_table=format_table(sanafe.framework_attributes),
    )

def generate_framework_page() -> str:
    return FRAMEWORK_PAGE_TEMPLATE.format(
        framework_table=format_table(sanafe.framework_attributes),
    )

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    framework_path = os.path.join(OUTPUT_DIR, "framework.rst")
    with open(framework_path, "w") as out_file:
        out_file.write(generate_framework_page())
    print(f"Wrote {framework_path}")

    for model_name in sanafe.model_attributes:
        out_path = os.path.join(OUTPUT_DIR, f"{model_name}.rst")
        with open(out_path, "w") as out_file:
            out_file.write(generate_model_page(model_name))
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
