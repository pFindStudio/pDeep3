import setuptools
from configparser import ConfigParser
from pkg_resources import parse_version
from sys import platform
assert parse_version(setuptools.__version__) >= parse_version('36.2')

config = ConfigParser(delimiters=['='])
config.read('settings.ini')
cfg = config['DEFAULT']

license_options = {
    'apache2': (
        'Apache Software License 2.0',
        'OSI Approved :: Apache Software License'
    ),
    'MIT': (
        'MIT License',
        'OSI Approved :: MIT License'
    )
}
status_options = {
    '1': 'Planning',
    '2': 'Pre-Alpha',
    '3': 'Alpha',
    '4': 'Beta',
    '5': 'Production/Stable',
    '6': 'Mature',
    '7': 'Inactive'
}
maximum_python3_available = 8

with open("requirements.txt") as requirements_file:
    requirements = []
    for line in requirements_file:
        line = line.strip()
        requirements.append(line)

setuptools.setup(
    name=cfg["lib_name"],
    license=license_options[cfg["license"]][0],
    classifiers=[
        f'Development Status :: {cfg["status"]} - {status_options[cfg["status"]]}',
        f'Intended Audience :: {cfg["audience"]}',
        f'License :: {license_options[cfg["license"]][1]}',
        f'Natural Language :: {cfg["language"]}',
    ] + [
        f'Programming Language :: Python :: 3.{i}' for i in range(
            int(cfg["min_python"].split(".")[1]),
            maximum_python3_available + 1
        )
    ],
    version=cfg["version"],
    description=cfg["description"],
    keywords=cfg["keywords"],
    author=cfg["author"],
    author_email=cfg["author_email"],
    url=cfg["url"],
    packages=setuptools.find_packages(),
    # TODO: Modifying this should allow to remove the MAINFEST.in
    include_package_data=True,
    install_requires=requirements,
    python_requires=f'>={cfg["min_python"]},<{cfg["max_python"]}',
    zip_safe=False,
)
