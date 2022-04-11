export GIT_VERSION=`git log -1 --pretty=format:"%H"`
if [[ $(git diff --stat) != "" ]];
	then GIT_VERSION="$GIT_VERSION-dirty"
fi

echo $GIT_VERSION;
