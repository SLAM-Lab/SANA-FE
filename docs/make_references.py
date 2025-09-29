#!/usr/bin/env python3
"""
Citation reference generator for SANA-FE project.
Reads references.bib and updates both README.md and docs/index.rst.

Usage (from docs/ directory):
    python make_references.py

Requirements:
    pip install bibtexparser
"""

import re
import sys
from pathlib import Path
try:
    import bibtexparser
    from bibtexparser.bparser import BibTexParser
except ImportError:
    print("Error: bibtexparser not found. Install with: pip install bibtexparser")
    sys.exit(1)


def load_bibtex(bib_file: Path) -> dict:
    """Load and parse BibTeX file."""
    if not bib_file.exists():
        print(f"Error: {bib_file} not found!")
        sys.exit(1)

    with open(bib_file, 'r', encoding='utf-8') as f:
        parser = BibTexParser(common_strings=True)
        bib_database = bibtexparser.load(f, parser=parser)

    return {entry['ID']: entry for entry in bib_database.entries}


def format_authors(authors_str: str) -> str:
    """Format authors string, handling 'and' separators."""
    authors = [author.strip() for author in authors_str.split(' and ')]
    if len(authors) == 1:
        return authors[0]
    elif len(authors) == 2:
        return f"{authors[0]} and {authors[1]}"
    else:
        return f"{', '.join(authors[:-1])}, and {authors[-1]}"


def generate_markdown_citation(entries: dict) -> str:
    """Generate markdown citation section for README.md."""
    # Use the main SANA-FE paper (boyle2025sanafe)
    main_entry = entries.get('boyle2025sanafe')
    if not main_entry:
        print("Warning: Main citation 'boyle2025sanafe' not found in BibTeX")
        return ""

    authors = format_authors(main_entry['author'])
    title = main_entry['title']
    journal = main_entry['journal']
    volume = main_entry['volume']
    number = main_entry['number']
    pages = main_entry['pages']
    year = main_entry['year']
    doi = main_entry['doi']

    page_start = pages.split("--")[0]
    page_end = pages.split("--")[1]
    citation_text = f"""# Citation

We hope that you find this project useful. If you use SANA-FE in your work,
please cite our paper:

{authors},\n"{title}," in\n{journal}, vol. {volume}, no. {number}, pp. {page_start}–{page_end}, {year},\n[doi:{doi}](https://doi.org/{doi}).

```bibtex
@article{{boyle2025sanafe,
  title={{{title}}},
  author={{{main_entry['author']}}},
  journal={{{journal}}},
  volume={{{volume}}},
  number={{{number}}},
  pages={{{pages}}},
  year={{{year}}},
  doi={{{doi}}}
}}
```"""
    return citation_text


def generate_markdown_references(entries: dict) -> str:
    """Generate markdown references section for README.md."""
    references = []

    for entry in entries.values():
        authors = format_authors(entry['author'])
        title = entry['title']
        year = entry['year']
        volume = entry.get('volume')
        number = entry.get('number')
        pages = entry.get('pages')
        year = entry['year']
        doi = entry['doi']

        if entry['ENTRYTYPE'] == 'article':
            journal = entry['journal']
            page_start = pages.split("--")[0]
            page_end = pages.split("--")[1]
            ref = f'{authors},\n"{title},"\nin {journal}, vol. {volume}, no. {number}, pp. {page_start}–{page_end}, {year},\n[doi:{doi}](https://doi.org/{doi}).'
        else:  # inproceedings
            booktitle = entry['booktitle']
            address = entry.get('address', '')
            if address:
                ref = f'{authors},\n"{title},"\nin {booktitle}, {address}, {year},\n[doi:{doi}](https://doi.org/{doi}).'
            else:
                ref = f'{authors},\n"{title},"\nin {booktitle}, {year},\n[doi:{doi}](https://doi.org/{doi}).'

        references.append(ref)

    references_text = "# References\n\n" + "\n\n".join(references)
    return references_text


def generate_rst_citation(entries: dict) -> str:
    """Generate RST citation section for index.rst."""
    main_entry = entries.get('boyle2025sanafe')
    if not main_entry:
        return ""

    title = main_entry['title']
    author = main_entry['author']
    journal = main_entry['journal']
    doi = main_entry['doi']
    volume = main_entry['volume']
    number = main_entry['number']
    pages = main_entry['pages']
    year = main_entry['year']
    doi = main_entry['doi']

    citation_text = f"""Citation
========

We hope that you find this project useful. If you use SANA-FE in your work,
please cite our paper:

.. code-block:: bibtex

   @article{{boyle2025sanafe,
     title={{{title}}},
     author={{{author}}},
     journal={{{journal}}},
     volume={{{volume}}},
     number={{{number}}},
     pages={{{pages}}},
     year={{{year}}},
     doi={{{doi}}}
   }}"""
    return citation_text


def generate_rst_references(entries: dict) -> str:
    """Generate RST references section for index.rst."""
    references = []

    for entry in entries.values():
        print(entry)
        authors = format_authors(entry['author'])
        title = entry['title']
        year = entry['year']
        doi = entry['doi']

        # Optional entries
        volume = entry.get('volume', None)
        number = entry.get('number', None)
        pages = entry.get('pages', None)

        if entry['ENTRYTYPE'] == 'article':
            journal = entry['journal']
            page_info = ""
            if pages:
                # Page(s) could be given as a range or single value
                if "--" in pages:
                    page_start = pages.split("--")[0]
                    page_end = pages.split("--")[1]
                    page_info = f"pp. {page_start}–{page_end}"
                else:
                    page_info = f"p. {pages.strip()}"

            # Build reference up with some optional fields
            ref_parts = [f'{authors},\n"{title},"']
            ref_parts.append(f'in {journal}')

            # Add volume and number if available
            if volume and number:
                ref_parts.append(f", vol. {volume}, no. {number}")
            elif volume:
                ref_parts.append(f", vol. {volume}")
            elif number:
                ref_parts.append(f", no. {number}")

            # Add page info if available
            if page_info:
                ref_parts.append(f', {page_info}')

            ref_parts.append(f', {year},\n`doi:{doi} <https://doi.org/{doi}>`_.')
            ref = ''.join(ref_parts)
        else:  # inproceedings
            booktitle = entry['booktitle']
            address = entry.get('address', '')
            if address:
                ref = f'{authors},\n"{title},"\nin {booktitle}, {address}, {year},\n`doi:{doi} <https://doi.org/{doi}>`_.'
            else:
                ref = f'{authors},\n"{title},"\nin {booktitle}, {year},\n`doi:{doi} <https://doi.org/{doi}>`_.'

        references.append(ref)

    references_text = "References\n==========\n\n" + "\n\n".join(references)
    return references_text


def update_file_section(file_path: Path, pattern: str, new_content: str, section_name: str):
    """Update a specific section in a file."""
    if not file_path.exists():
        print(f"Warning: {file_path} not found, skipping.")
        return False

    content = file_path.read_text(encoding='utf-8')

    if not re.search(pattern, content, re.DOTALL):
        print(f"Warning: {section_name} section not found in {file_path}")
        return False

    updated_content = re.sub(pattern, new_content + '\n\n', content, flags=re.DOTALL)
    file_path.write_text(updated_content, encoding='utf-8')
    print(f"Updated {section_name} in {file_path}")
    return True


def main():
    """Main function to update all citation files."""
    # File paths (relative to docs/ directory)
    bib_file = Path("../references.bib")
    readme_file = Path("../README.md")
    rst_file = Path("index.rst")

    print("Loading citations from references.bib...")
    entries = load_bibtex(bib_file)
    print(f"Found {len(entries)} citations")

    # Update README.md
    if readme_file.exists():
        print("\nUpdating README.md...")

        # Update citation section
        citation_pattern = r'# Citation.*?(?=# [A-Z]|\Z)'
        new_citation = generate_markdown_citation(entries)
        update_file_section(readme_file, citation_pattern, new_citation, "Citation")

        # Update references section
        references_pattern = r'# References.*?(?=# [A-Z]|\Z)'
        new_references = generate_markdown_references(entries)
        update_file_section(readme_file, references_pattern, new_references, "References")

    # Update index.rst
    if rst_file.exists():
        print("\nUpdating docs/index.rst...")

        # Update citation section
        citation_pattern = r'Citation\n========.*?(?=\n[A-Z][a-z]*\n=+|\Z)'
        new_citation = generate_rst_citation(entries)
        update_file_section(rst_file, citation_pattern, new_citation, "Citation")

        # Update references section
        references_pattern = r'References\n==========.*?(?=\n[A-Z][a-z]*\n=+|\Z)'
        new_references = generate_rst_references(entries)
        print(new_references)
        update_file_section(rst_file, references_pattern, new_references, "References")

    print("\nCitation update complete!")
    print("\nNext steps:")
    print("1. Review the updated files")
    print("2. Commit changes to git")
    print("3. To add new citations, edit references.bib and rerun this script")

if __name__ == "__main__":
    main()