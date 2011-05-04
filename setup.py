from setuptools import setup, find_packages

version = '0.1dev'

setup(name='fragions',
      version=version,
      description="Scripts to calculate frequency of fragment ions",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='Danny Navarro',
      author_email='j@dannynavarro.net',
      url='http://dannynavarro.net',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=['numpy']
          # -*- Extra requirements: -*-
      ],
      entry_points="""\
        [console_scripts]
          exp_db.py = fragions.exp_db:main
          dump_csv.py = fragions.dump_csv:main
      """
      )
