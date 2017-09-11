.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/Biosustain/cobra_gecko/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

cobra gecko could always use more documentation, whether as part of the
official cobra gecko docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/Biosustain/cobra_gecko/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `cobra_gecko` for local development.

1. Fork the `cobra_gecko` repo on GitHub.
2. Clone your fork locally::

    git clone git@github.com:your_name_here/cobra_gecko.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    mkvirtualenv cobra_gecko
    cd cobra_gecko/
    python setup.py develop

4. Create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox::

    tox

   Or do so in parallel with ``detox``. To get tox and detox, just pip install them into your virtualenv.

6. Commit your changes using `semantic commit messages <https://seesparkbox.com/foundry/semantic_commit_messages>`__ and push your branch to GitHub::

    git add .
    git commit -m "feat: your detailed description of your changes"
    git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.
