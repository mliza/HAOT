from setuptools import setup, find_packages

setup(
    name="hoat",
    version="1.0.0",
    author="Martin E. Liza",
    author_email="mliza1191@gmail.com",
    description="This package contains tools used for aero-optics post
    processing",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mliza/HAOT",
    packages=find_packages(),
    license="MIT",
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.6.0",
        "molmass>=2024.10.25"
    ],
    entry_points={
        "console_scripts": [
            "hoat-cli=hoat.__main__:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
