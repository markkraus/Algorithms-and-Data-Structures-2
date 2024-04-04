# Repository Guide
Link: https://drive.google.com/file/d/1Mf0Rp1QxbXa_pkKKRHfb9n28WYlffDxW/view?usp=sharing<br>

## Project Documentation
https://github.com/cs1501-2241/cs1501-assignment3-markkraus<br>

# LZW.java
- Array & Trie Codebooks<br>
- Standard input/output reading/writing<br>
- Dynamic dictionary resetting<br>

# How to Run
javac *.java<br>
To Compress without a dictionary reset: java LZW - n < _test file.type_ > _test file.type_.lzw <br>
To compress with a dictionary reset:    java LZW - r < _test file.type_ > _test file.type_.lzw <br>
To expand the compressed file:          java LZW + < _test file.type_.lzw > _test file.type_.rec <br>
