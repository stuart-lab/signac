# Contributing to Signac development

The goal of this guide is to help you get up and contributing to Signac as 
quickly as possible. The guide is divided into two main pieces:

1. Filing a bug report or feature request in an issue.
1. Suggesting a change via a pull request.

Please note that Signac is released with a [Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this project, 
you agree to abide by its terms.

## Issues

GitHub issues are intended for filing bug reports, feature requests, or
requests for changes to the documentation. For general discussion of analysis
steps, opening a discussion topic on the Signac
[discussion forum](https://github.com/stuart-lab/signac/discussions) is preferred.

When filing a bug report, it is extremely useful if you can include a minimal 
reproducible example so that we can quickly verify the problem, and then figure 
out how to fix it. There are three things you need to include to make your 
example reproducible: packages, data, code.

1.  Include the output of `sessionInfo()` in your issue. This lets us see what
    **packages** you have loaded and their version number
  
1.  The easiest way to include **data** is to use the small
    example dataset included in Signac (`atac_small`). Sometimes bugs are more
    complex and can't be reproduced using the small example data. In these
    cases, try to use one of the datasets used in the Signac vignettes. These 
    are all publicly available data, and instructions on how to download these
    datasets are available on the vignette pages.
  
1.  Spend some time ensuring that your **code** is easy for others to read:
    
    * do your best to remove everything that is not related to the problem.  
     The shorter your code is, the easier it is to understand.
  
    * make sure you've used spaces and your variable names are concise, but
      informative
  
    * use comments to indicate where your problem lies
    
    * properly format the code in your issue. Information about formatting is
    available [here](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax#quoting-code)
    
Before posting your issue, make sure to start a new R session and run your
reproducible example, and verify that it reproduces your issue.

## Pull requests

To contribute a change to Signac, you follow these steps:

1. Fork the Signac repository
1. Clone your forked version of the repository
1. Run `git checkout develop` to change to the develop branch
1. Create a new branch for your changes (`git checkout -b <new branch name>`)
1. Add your changes to the branch (`git commit <files to commit`)
1. Push your branch to github (`git push`)
1. Open a pull request (PR) in the Signac repository

If you get stuck at any point, please reach out for help by opening a discussion
topic on the Signac GitHub page.

If you're not familiar with git or GitHub, you can find some useful information
[here](https://github.com/git-guides).

When reviewing your PR, we will check for the following things:

1.  __Motivation__. Your pull request should clearly and concisely motivate the
    need for change. You need to describe the problem and show
    how your pull request solves it. If the PR fixes a bug or implements a 
    feature request that was reported in an issue, you should reference that
    issue in your PR. 

1.  __Only related changes__. Before you submit your pull request, please
    check to make sure that you haven't accidentally included any unrelated
    changes. These make it harder to see exactly what's changed, and to
    evaluate any unexpected side effects.

1.  If you're adding new parameters or a new function or updating function 
    __documentation__, you'll also need to update the documentation pages by
    running [roxygen](https://github.com/klutometis/roxygen).
    Make sure to re-run `devtools::document()` on the code before submitting.

1.  We encourage you to add a [testthat](https://github.com/r-lib/testthat)
    __unit test__ if including new functionality.

This can seem like a lot of work but don't worry if your pull request isn't
perfect. We are very happy to receive contributions and will help to update your
PR to make it acceptable to merge into the Signac package. If you can't include
all the things listed above, please still submit your PR and we can work
together on improving it.

Finally, remember that Signac is being used by many people, and so we try not to
add changes that will alter the existing functionality unless absolutely
necessary, as that may break someone's code. Cases where existing functionality
will be changed will be reserved for major version releases.
