VERSION=`grep "define version" src/*cpp | awk '{print substr($3,2,length($3)-2)}'`
echo BUILDING ${VERSION}
cd src
make clean
make release
cd ../
mkdir -p duohmm_${VERSION}/
cp src/duohmm duohmm_${VERSION}/
cp python/mapavg.py  duohmm_${VERSION}/
cp -r example duohmm_${VERSION}/
tar -pczf duohmm_${VERSION}.tar.gz duohmm_${VERSION}
