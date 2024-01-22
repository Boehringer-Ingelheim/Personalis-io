# Contribution guidelines

After making changes, you'll need to take the following steps:

-   Rebuild roxygen: `Rscript -e "devtools::document()"`
-   Run R CMD CHECK: `Rscript -e 'devtools::check()'`
-   Build documentation website: `Rscript  -e 'pkgdown::build_site_github_pages(install=TRUE)'`
-   Push the `/docs` folder to the `gh-pages` branch of the repo to update the documentation website

## Making a release

-   bump the version according to semver.org
-   make a release on GitHub

After making the release, bump the version to x.x.x.9000 to mark it as a development version
