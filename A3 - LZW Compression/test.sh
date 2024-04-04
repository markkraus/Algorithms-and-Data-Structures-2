#! /bin/sh
# Usage: ./test.sh <file to compress>
# Use any of the files provided inside the "TestFiles" folder
# Example: ./test.sh TestFiles\code.txt
echo "javac *.java"
javac *.java
echo "java LZW - < $1 > $1.lzw"
java LZW - < $1 > $1.lzw
echo "java LZW + < $1.lzw > $1.rec"
java LZW + < $1.lzw > $1.rec
echo "diff $1 $1.rec"
diff $1 $1.rec