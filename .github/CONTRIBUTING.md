## Contributor guidelines

First of all, thank you for contributing to GECKO!! Anybody is welcome to contribute, but please abide by the following guidelines.

You can contribute in 2 main ways: by creating issues, and by sending pull requests (PRs) with corrections or additions to the source code.

### Reporting issues

[Report an issue](https://github.com/SysBioChalmers/GECKO/issues) if you note any of the following:

* Missing feature you would like GECKO to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

If you are unsure about the issue, consider asking first in our [Gitter chat room](https://gitter.im/SysBioChalmers/GECKO).

When creating the issue, please make sure:

* You tested your code (if any) with all of GECKO's [requirements](https://github.com/SysBioChalmers/GECKO).
* You did your analysis in the `main` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.
* You checked that a similar issue does not exist already

Feel free to also comment on any of the [open issues](https://github.com/SysBioChalmers/GECKO/issues). When doing so, please comply with our [code of conduct](https://github.com/SysBioChalmers/GECKO/blob/main/.github/CODE_OF_CONDUCT.md).

Finally, if you like GECKO please remember to 'star' our Github page (click on the star at top right corner), that way we also have an idea of who is using GECKO!

### Contributing to GECKO

Do you want to contribute to GECKO with some additions or improvements? Consider starting by raising an issue and assign it to yourself to describe what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could start with. Also, feel free to browse our [open issues](https://github.com/SysBioChalmers/GECKO/issues): Anything tagged with "help wanted" is open to whoever wants to implement it!

Here's how to set up GECKO for local development to contribute smaller features or changes that you can implement yourself:

1. First of all, make sure that you have all [requirements](https://github.com/SysBioChalmers/GECKO) for running GECKO.

2. [Fork the GECKO repository](https://github.com/SysBioChalmers/GECKO/fork) on GitHub.

3. Clone your fork locally:
    ```
    git clone https://github.com/<your Github name>/GECKO.git
    ```

4. Check out the branch that you want to contribute to. Most likely that will be ``develop``:
    ```
    git checkout develop
    ```

5. Create a branch for local development based on the previously checked out branch ([see below](#branching-model) for details on the branching model and how to name your branch):
    ```
    git checkout -b name-of-your-branch
    ```

6. Now you can make your changes locally!
  * New scripts (if any) should start with a commented section describing the script and explaining the inputs/outputs. If you are uncertain on the style to follow, check out pre-existing functions.
  * Try to document as much as possible; this will help the review process.

7. Commit your changes and push your branch to GitHub.
    ```
    git add .
    git commit -m "Title of your commit"
    git push origin name-of-your-branch
    ```
    [See below](#semantic-commits) for recommendations on how to name your commits. In case of larger updates, you can of course make several commits on a single contribution. However, if you need to do too many commits, consider if your contribution could be instead split into separate PRs (making it easier for reviewing later).

8. Submit a PR through the GitHub website (https://help.github.com/articles/creating-a-pull-request-from-a-fork/) to the `develop` branch of the original SysBioChalmers repo (not your fork). We recommend ticking the box "Allow edits from maintainers" if you wish for us to be able to contribute directly to your branch (speeding-up the reviewing process).

Note that steps 3, 4, 5 and 7 can be done, if you prefer, with any git client, such as [Github Desktop](https://desktop.github.com/).

Try to keep pull requests relatively small: they should ideally provide one new feature or fix one bug. A few small bug-fixes maybe combined, but avoid packing pull requests with too many diverse and somewhat unrelated development. Smaller pull requests reduce reviewing time and make sure that code updates will make it into new releases faster.

Finally, and for larger features that you want to work on collaboratively, you may consider to first request to join our development team to get write access to the repository so that you can create a branch directly in the main repository. On this new branch, you can push your changes directly to the main repository and when finished, submit a pull request from that branch to `develop`. [See below](#development-team-guidelines) for more details.

Thank you very much for contributing to GECKO!

#### Branching model

* `develop`: Is the branch all pull-requests should be based on.

* `main`: Is only touched by the administrator and is the branch with the tested & reviewed codebase. Each merge commit to it is associated to a release.

* `{chore, doc, feat, fix, refactor, style}/descriptive-name`: Any other branch created in the model. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/bug` or `feat/new-algorithm`. [See below](#semantic-commits) for more details on the possible actions you can use.

#### Semantic commits

Please use concise descriptive commit messages. Ideally, use [semantic commit messages](http://karma-runner.github.io/2.0/dev/git-commit-msg.html) to make it easier to show what you are aiming to do:

`action: brief description`

`action` refers to what exactly are you doing in the commit:
* `chore`: updating toolbox, data files, etc.
* `doc`: updating documentation or explanatory comments in functions.
* `feat`: new feature added.
* `fix`: something that was incorrect and now has been corrected.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes to functions or data.

More examples [here](https://github.com/SysBioChalmers/GECKO/commits/master). A more detailed explanation or comments is encouraged to be left in the commit description.

## Development team guidelines

This section is meant for the development team of GECKO. As a member of the development team, you should comply with all [previous contributor guidelines](#contributor-guidelines) as well. Besides that, please also consider the following guidelines.

### Creating pull requests

Changes to the source code _may not_ be directly committed to the `main` or `develop` branches (in fact they are protected). Commits are made to side-branches, after which pull requests are made for merging with `develop`. For this, follow the [instructions](#contributing-to-the-model) for contributors, but consider that as members of the development team have write access to the repository, you can create a branch directly in the main repository without needing to fork, for your convenience. This means that you can:

* Skip step 2 of the contribution process.
* In step 3 of the contribution process, clone directly the original repo:
  ```
  git clone https://github.com/SysBioChalmers/GECKO.git
  ```

Follow all other steps in the same way. Also, when creating your pull request (or after you have created it):

* Choose 1 or more members of the team (ideally with knowledge on the pull request) as reviewers. Note that the person making the pull request and the reviewer _cannot_ be the same person.
* Assign appropriate [labels](https://github.com/SysBioChalmers/GECKO/issues/labels).

### Reviewing pull requests

Every pull request must be approved by at least one reviewer before it can be merged. When reviewing someone else's pull request, keep in mind the following aspects:

* **Functionality:** The automated tests that are run as GitHub Actions should pass without any problems. Have you written a new function? Make sure that it is also tested in `tests/unit_tests/geckoCoreFunctionTests.m`.
* **Documentation:** The reasoning for each modification should be provided, as documentation or as a comment in the pull request.
* **Reproducibility:** If there are any added scripts, make sure that if you run them, they achieve the desired purpose.
* **Style:** Ensure that changes to the codebase have a compliant style, and new datasets (if any) are straight-forward to understand.
* When commenting in the review, please comply with our [code of conduct](https://github.com/SysBioChalmers/GECKO/blob/main/.github/CODE_OF_CONDUCT.md).
* Avoid vague comments and try to be as explicit as possible (e.g.: _"Please include X here"_ instead of _"X could be included here"_).
* As much as possible, try to keep the review process in the pull request discussion, and not in separate private emails.

## Administrator guidelines

This section is meant for the administrator of this repo. The main duties of the administrator are:

* To keep the repository clean and organized, i.e. avoid redundancy in functions and/or data, and keep coherency in naming of files.
* To help in the reviewing process of external pull requests by assigning reviewers and [labels](https://github.com/SysBioChalmers/GECKO/issues/labels), if applicable.
* To keep [issues](https://github.com/SysBioChalmers/GECKO/issues) with the proper labels, and to close issues once they are fixed in the `main` branch.
* In cases of disagreement between contributors, to decide how to resolve the issue.
* To merge open pull requests into `develop` regularly (see [below](#merging-contributions)).
* To generate new releases of the model regularly (see [below](#releasing-a-new-version)).

### Merging contributions

The following points should be considered when merging branches to `develop`:

* Make sure the branch gets accepted by at least one developer with writing access.
* Wait at least a day before merging, to allow other developers to inspect the pull request.
* To simplify the git history, use "Squash and merge" in pull requests to `develop` (but **never** use this option when merging to `main`, see also below).
* As soon as the branch is merged, confirm that `develop` is still possible to merge to `main` (this can be checked [here](https://github.com/SysBioChalmers/GECKO/compare/develop)). If conflicts appear, fix the conflict _locally_ as soon as possible in `develop` and then push it (note, **DO NOT** pull any other changes from `main` to `develop`, just the single file that is creating the conflict).

### Releasing a new version

* A merge of `develop` with `main` invokes a new release.
* A new release should be made as soon as there is substantial new work in `develop` (as rule of thumb, after 3-4 pull request merges or ~20 commits).

GECKO follows [semantic versioning](https://semver.org/):

* A `major` release is when substantial new features are added and backwards compatibility is broken.
* A `minor` release is when substantial new features are added and backwards compatibility is not broken.
* A `patch` release is the most common one and is done when only few things have changed, such as:
  - Small new features.
  - Bug fixes.
  - Updating toolboxes.
  - Re-organization of data
  - Refactoring of code.

When releasing, please follow these steps:

1. Create a pull request from `develop` to `main`:
   - Specify the intended version in the title, e.g. `GECKO 3.0.1`
   - Indicating all new features/fixes/etc. and referencing every previous pull request included (examples [here](https://github.com/SysBioChalmers/GECKO/releases)).
   - If any [issue](https://github.com/SysBioChalmers/GECKO/issues) gets solved in the release, write in the pull request description "Closes #X", where "X" is the issue number. That way the issue will be automatically closed after merge.
2. Wait at least a day for at least one approval. The GitHub Actions must also pass successfully.
3. Merge the pull request from `develop` to `main` with the option "Create a merge commit", **not** "Squash and merge".
4. Edit the `version.txt` in the `main` branch directly [link](https://github.com/SysBioChalmers/GECKO/edit/main/version.txt) to refer to the new version number. Commit the change directly to the `main` branch.
5. Make the [new release](https://github.com/SysBioChalmers/GECKO/releases/new) at GitHub:
   - Define a new tag with the version, prefixed by `v`: e.g. `v3.0.1`.
   - Copy the description from the corresponding PR from step 1.
