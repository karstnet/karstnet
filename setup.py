import setuptools

# # Set long_description
# with open("README.md", "r") as file_handle:
#     long_description = file_handle.read()

# Load version
__version__ = '0.0.0' # default
with open(f'src/karstnet/_version.py', 'r') as f:
    exec(f.read())

# ----------------------------------------------------------------------------
setuptools.setup(
    version=__version__,
)
