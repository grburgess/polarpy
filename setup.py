import os

from setuptools import setup
import versioneer

# Create list of data files
def find_data_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join('..', path, filename))

    return paths


extra_files = find_data_files('polarpy/data')


setup(
    version=versioneer.get_version()
     cmdclass=versioneer.get_cmdclass()

    include_package_data=True,
)
