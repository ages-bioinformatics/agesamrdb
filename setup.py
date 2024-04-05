from setuptools import setup


setup (
      name = "agesamrdb",
      version = 1.1,
      description = "",
      author = "Patrick Hyden, Tobias Mösenbacher",
      author_email = "patrick-christian.hyden@ages.at",
      packages = ["agesamrdb"],
      data_files = [("agesamrdb", ["agesamrdb/data/phenotypes_classnames.tsv"])],
      include_package_data=True,
      install_requires=["pandas", "sqlalchemy>=2.0", "mysql-connector-python"],
      scripts = ["amrdb_add_results.py", "update_resfinder_database.py"],
      long_description = "This tool allows to create a relational database" \
              + " from ResFinder and store ResFinder (PointFinder) results" \
              + " alongside with additional tools results: currently " \
              + "implemented: Bakta, ISEScan2, Plasmidfinder, AMRFinderPlus",
      license = "MIT",
      platforms = "Linux, Mac OS X"
)

