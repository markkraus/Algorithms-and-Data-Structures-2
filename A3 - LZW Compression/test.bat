:: Usage: test.bat <file to compress>
:: example: test.bat TestFiles\code.txt
:: Use any of the files provided inside the "TestFiles" folder
@echo off
javac *.java
java LZW - < %1 > %1.lzw
java LZW + < %1.lzw > %1.rec
fc %1 %1.rec