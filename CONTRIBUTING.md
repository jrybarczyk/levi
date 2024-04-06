
# Contributing to levi

We are thrilled to have you contribute to our project! Our goal is to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing or implementing new features
- Interested in becoming a maintainer

## How to Report Bugs, Seek Help, or Engage in Discussions

As levi is under development, we currently do not have a dedicated forum for help. Please submit your questions, bugs, and discussion topics to our [GitHub issues page](https://github.com/jrybarczyk/levi/issues).

For bug reports, please include as much of the following information as possible to help us address the issue quickly:

- A concise summary and/or background
- Steps to reproduce the issue. Be specific! **Include sample code if possible.**
- Your expected outcome vs. what actually happened
- Any error messages or stack traces you encountered
- Notes or hypotheses on why the issue might be occurring, or any attempted solutions that didn't work

Thorough bug reports are greatly appreciated as they allow our development team to quickly identify and resolve issues. Once we verify a bug, it will be added to our GitHub issues for tracking and resolution.

## How to Contribute Code Changes via GitHub

We utilize GitHub for code hosting, issue tracking, and accepting pull requests.

Pull requests are the best way to propose changes to the codebase. We recommend reading through the [GitHub flow](https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow) to understand this process better.

The basic steps for creating a pull request are:

- Fork the repository and create your branch from master.
- Make your changes, commit them to your branch, and then push them to your GitHub fork.
- Once you're done with your changes, navigate to your fork on GitHub and initiate a Pull Request. It will automatically update if you need to make further adjustments.

## Guidelines for a Successful Pull Request

Here are some tips to help make your pull requests successful:

- Adopt the [Bioconductor style guide](https://contributions.bioconductor.org/develop-overview.html) for your R code. This helps keep the codebase consistent and easy to read.
- Use spaces (4) for indentation rather than tabs.
- If necessary, update the documentation to reflect your changes.
- **Write tests** for new features. Quality tests are crucial for maintaining robust code. We encourage using the `testthat` package in R for testing. Look at the existing tests for inspiration or refer to the [R packages book](https://r-pkgs.org/tests.html) for guidance on writing good tests.
- Be aware that your contributions will be licensed under the same terms as the project.

When you submit your pull request, automated tests will run as part of our continuous integration setup. We appreciate constructive discussion on the best practices for coding and welcome comments on your pull request as a forum for such exchanges.
