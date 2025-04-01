import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="seqranker",
    version="0.1.4",
    author="Till-Hendrik Macher",
    author_email="macher@uni-trier.de",
    description="dbDNA - A phylogeny- and expert identifier-driven grading system for reliable taxonomic annotation of (meta)barcoding data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TillMacher/dbdna",
    packages=setuptools.find_packages(),
    license="MIT",
    install_requires=[
        'Bio>=1.7.1',
        'biopython>=1.84',
        'ete3>=3.1.3',
        'joblib>=1.4.2',
        'numpy>=2.2.4',
        'pandas>=2.2.3',
        'playwright>=1.51.0',
        'plotly>=5.9.0',
        'PyPDF2>=3.0.1',
        'Requests>=2.32.3',
        'requests_html>=0.10.0',
        'tqdm>=4.66.4',
        'xmltodict>=0.13.0',
        'lxml_html_clean>=0.4.1',
        'openpyxl>=3.1.5',
        'pyarrow>=19.0.1',
        'fastparquet>=2024.11.0',
        'kaleido>=0.2.1',
    ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    entry_points={
        "console_scripts": [
            "seqranker = seqranker.__main__:main",
        ]
    },
)
