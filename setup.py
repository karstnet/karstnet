import setuptools

with open("README.md", "r") as file_handle:
    long_description = file_handle.read()

setuptools.setup(
    name='karstnet',
    version='1.2.1',
    author="Philippe Renard & Pauline Collon",
    description="Statistical characterization of karstic network",
    long_description=long_description,
    packages=setuptools.find_packages()
)
