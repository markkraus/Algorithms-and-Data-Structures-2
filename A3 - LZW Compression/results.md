# Compression Table
| Compressed File Name | Original Uncompressed File Size | 12-bit Codewords | Variable-bit Codewords w/ Dictionary Reset | Variable-bit Codewords w/out Dictionary Reset | Predefined Unix Compress Program |
| :--------------------------: | :---------------------------------------------: | :------------------------------------------------------------------: | :----------------------------: | :--------------------------------: | :-----------: |
|   all.lzw              |3,031,040 bytes | 1,846,855 bytes              | **24,290 bytes**  |1,792,781 bytes| 1,179,467 bytes              | 
|   assig2.lzw           | 87,040 bytes   | 74,575 bytes              | **40,039 bytes**  | **40,039 bytes**  | 40,040 bytes              |
|   Lego-big.lzw         | 93,371 bytes   |128,974 bytes               | **122,493 bytes** |**122,493 bytes**| n/a              |
|   bmps.lzw             |1,105,920 bytes | 925,079 bytes              | **80,913 bytes**  | **80,913 bytes**  | **80,913 bytes**              |
|   code2.lzw            | 55,173 bytes   | 23,776 bytes              | **20,318 bytes**  | **20,318 bytes**  | 20,319 bytes              |
|   code.lzw             | 69,516 bytes   |30,853 bytes                | **24,290 bytes**  | **24,290 bytes**  | 24,291 bytes              |
|   edit.lzw             | 236,328 bytes  | 250,742 bytes              | 152,230 bytes | 156,409 bytes | **151,111 bytes**               |
|   frosty.lzw           | 126,748 bytes  |  177,454 bytes             | 171,169 bytes|**163,789 bytes**| n/a              |
|   gone_fishing.lzw     |  17,336 bytes  | 9,278 bytes              | **8,962 bytes**   |   **8,962 bytes** | 8,964 bytes               |
|   large.lzw            |1,202,961 bytes | 599,407 bytes              | 523,714 bytes | **497,971 bytes** |519,465 bytes               |
|   medium.lzw           |  24,599 bytes  | 13,076 bytes              | **12,439 bytes**  | **12,439 bytes**  | 12,440 bytes              |
|   texts.lzw            |1,382,400 bytes |1,012,180 bytes               | 590,558 bytes | 597,847 bytes |  **589,697 bytes**             |
|   wacky.lzw            | 921,654 bytes  | 4,303 bytes              | **3,951 bytes**   | **3,951 bytes**   |  3,952 bytes             |
|   winnt256.lzw         | 157,044 bytes  |  159,050 bytes   | **62,931 bytes**  | **62,931 bytes**  |  **62,931 bytes**             |

## Reflection
After compressing all of the test files with the 4 variations of compressing, it is clear that the variable-bit codewords gets the job done most efficiently. What's also interesting to note is that the two different variations of 
variable-bit codewords implementations usually compressed to the same amount of bytes. I speculate this happens because not all of the files will hit the codeword max limit, and so the files in which the two variations 
have equal amount of compressed bytes, I believe there was never a need to decide if the codebook needed to be reset. The files where there was a difference, however, showed how resetting the dictionary could change the
compression size.

## Best Ratios
The best compression ratio for the 12-bit Codewords was **wacky.lzw (4,303:921,654)**.<br>
- I think for a .bmp file, the tight bound on codebook size may help because the patterns gotten in this file type will be minimized and repeated.<br>
<br>
The best compression ratio for the Variable-bit Codewords w/ Dictionary Reset was **all.lzw (24,290:3,031,040)**.<br>
- Because the file size is so large, resetting can help introduce more new patterns and better compress the file over time.<br>
<br>
The best compression ratio for the Variable-bit Codewords w/out Dictionary Reset was **bmps.lzw (80,913:1,105,920)**.<br>
- I think with this file type, not resetting may help because the patterns may remain the same throughout the entire file.<br>
<br>
The best compression ratio for the Predefined Unix Compress Program was **wacky.lzw (3,952:921,654)**.<br>
- I think for the same reasons as the 12-bit codewords algorithm apply for this algorithm.<br>

## Worst Ratios
The worst compression ratio for the 12-bit Codewords was **Lego-big.lzw**.<br>
- I assume because a .gif file has many random patterns, the compression is hard to get with many reusable patterns, so the max codebook length is hit early.<br>
<br>
The worst compression ratio for the Variable-bit Codewords w/ Dictionary Reset was **Lego-big.lzw**.<br>
- Again, my reasoning is the same as for the algorithm above. I also think that even with a reset of the dictionary we don't find alot of patterns.<br>
<br>
The worst compression ratio for the Variable-bit Codewords w/out Dictionary Reset was **Lego-big.lzw**.<br>
- My reasoning is consistent with the reason I gave for the 12-bit Codeword algorithm.<br>
<br>
The worst compression ratio for the Predefined Unix Compress Program was **Lego-big.lzw**.<br>
- The compression was not even able to take place in this test.<br>

# In closing
Seeing the compression take place in different algorithms was interesting, because I can now see why it is difficult to determine the correct compression algorithm for a given file. Each file is incredibly unique and many factors play in to the compression of a file, like our codebook size, or how we want to handle a full codebook, and how to properly adjust the codebook as to not throw the program out of sync. This project was certainly a great way for me to enhance my knowledge on compression as it pertains to many different file types.
