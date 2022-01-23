import pathlib
from setuptools import setup

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(
    name="pcr_strainer",
    version="0.2.4",
    description="A tool for assessing the inclusivity of primer and probe oligonucleotides from diagnostic qPCR assays and amplicon sequencing schemes.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/KevinKuchinski/",
    author="Kevin Kuchinski",
    author_email="kevin.kuchinski@bccdc.ca",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["pcr_strainer"],
    include_package_data=True,
    install_requires=[
        'pandas',
    ],
    entry_points={
        "console_scripts": [
            "pcr_strainer=pcr_strainer.pcr_strainer:main",
        ]
    },
)
