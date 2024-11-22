#!/bin/sh
export TESTDIR=/var/lib/jenkins/workspace/PR_DOCS_EN
export LANG=en_US.UTF-8
make html
make latexpdf
rsync -cvl --delete -r ./build/html/ /var/www/html/docs/en
cp build/latex/SALMONdocument.pdf /var/www/html/docs/en
cp -r ./build/html/ .
cp build/latex/SALMONdocument.pdf .
git checkout -b pages
git add -A
git commit -m "modified by Jenkins" 
git push -f git@github.com:SALMON-TDDFT/SALMON-DOCS.git pages