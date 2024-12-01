from setuptools import setup, find_packages
import os

# Function to read the README file (for long description)
def read_readme():
    with open("README.md", "r") as f:
        return f.read()

# Function to gather all documentation files
def get_docs():
    docs_files = []
    for root, _, files in os.walk("docs"):
        for file in files:
            docs_files.append(os.path.relpath(os.path.join(root, file), start="docs"))
    return docs_files


setup(
    name="haot",
    version="1.0.0-b5",
    author="Martin E. Liza",
    author_email="mliza1191@gmail.com",
    description="Hypersonic Aerodynamic Optics Tools",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mliza/HAOT",
    packages=find_packages(),
    license="MIT",
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.6.0",
        "molmass>=2024.10.25",
        "ambiance>=1.3.1"
    ],
    entry_points={
        "console_scripts": [
            "hoat-cli=hoat.__main__:main",
        ],
    },
    include_package_data=True,
    package_data={
        "haot": get_docs(),  # Include documentation files in the package
    },
    zip_safe=False,
)
