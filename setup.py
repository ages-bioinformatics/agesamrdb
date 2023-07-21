from setuptools import setup


setup (
      name = "amrdb",
      version = 0.1,
      description = "",
      author = "Patrick Hyden",
      author_email = "patrick.hyden@ages.at",
      packages = ["amrdb"],
      package_dir = {"amr": "amr"},
      install_requires=["pandas", "sqlalchemy>=2.0", "mysql-connector-python"],
      scripts = ["amrdb_add_results.py", "update_resfinder_database.py"],
      long_description = """This tool allows to create a relational database from ResFinder and store ResFinder (PointFinder) results
      alongside with additional tools results: currently implemented: Bakta and ISEScan2""",
      license = "",
      platforms = "Linux, Mac OS X"
)

