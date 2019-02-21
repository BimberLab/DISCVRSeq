#!/bin/bash

# Adapted from: https://benlimmer.com/2013/12/26/automatically-publish-javadoc-to-gh-pages-with-travis-ci/

if [ "$TRAVIS_REPO_SLUG" != "BimberLab/DISCVRSeq" ];then
    echo -e "Incorrect fork, skipping doc publish: $TRAVIS_REPO_SLUG"
    exit 0
fi

if [ "$TRAVIS_PULL_REQUEST" != "false" ];then
    echo -e "Skipping publish for pull request"
    exit 0
fi

if [ "$TRAVIS_BRANCH" != "master" ];then
    echo -e "Skipping publish for branch $TRAVIS_BRANCH"
    exit 0
fi

echo -e "Publishing docs..."

git config --global user.email "travis@travis-ci.org"
git config --global user.name "travis-ci"
git clone --quiet --branch=gh-pages https://${GH_TOKEN}@github.com/BimberLab/DISCVRSeq gh-pages > /dev/null

# Commit and Push the Changes
cd gh-pages
git rm -rf ./*
cp -Rf ../build/docs/htmlDoc/* ./
git add -f .
git commit -m "Docs on travis build $TRAVIS_BUILD_NUMBER auto-pushed to gh-pages"
git push -fq origin gh-pages > /dev/null

echo -e "Published docs to gh-pages."