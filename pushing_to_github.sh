# For commiting to github for the first time
rm -rf .git
git init
# upload to github
git add .
git commit -m 'updated repository for viewORCA v0.2.0'
# For commiting to github for the first time
git branch -M main
git remote add origin git@github.com:geoffreyweal/viewORCA.git
git push -uf origin main
# push your new commit:
git push -u origin main
