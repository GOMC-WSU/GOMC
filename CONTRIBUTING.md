To contribute to GOMC please follow this procedure:

1)  Fork the GOMC repository

2)  Clone your forked repository

        git clone https://github.com/yourusername/GOMC.git

2)  Checkout the development branch:

        git checkout development

3)  Create a new branch based of the development branch named after your new feature

        git checkout -b nameofyourbranch

4)  Development on this branch checking for consistent results with the development branch by using the GOMC_Examples.py script

5)  When you are ready to add to GOMC, rebase your branch with the development branch.

        git remote add upstream https://github.com/GOMC-WSU/GOMC.git
        git rebase upstream/development

6)  Retest with the GOMC_Examples.py script

7)  If all the tests passed, open a pull request and select reviewers.
