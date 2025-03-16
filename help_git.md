0. Create .gitignore file and add ignored files
1. Create repository on github, nut no readme.md!
2. Run commands in project directory in terminal:

~~~bash
 git init
 git add . # (add all files)
 git commit -m "first commit"
 git branch -M main # rename master to main because of DEI :)
 git remote add origin git@github.com:osmianskiene/snakemake_2.git # Different for every project
 git push -u origin main
~~~

## Work

Commit & Push

or 

Commit -> Sync Changes

**If VSC hangs**, run from terminal:

~~~bash
git status
git add.
git commit -m "commit name"
~~~

If **Permission denied** when running star i.e.

~~~bash
chmod +x star/STAR
~~~