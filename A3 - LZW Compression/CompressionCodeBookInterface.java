/**
 * An interface for a codebook to be used for LZW compression */

public interface CompressionCodeBookInterface {

  /**
   * 
   * @return the current value of the codeword width
   */
  public int getCodewordWidth();

  /**
   * Adds a new string to the codebook in O(alphabet size*str.length()) time
   * @param str the String to be added to the codebook
   * @param flushIfFull a boolean indicating whether the codebook should be 
   *                    flushed when it is at its maximum capacity
   */
    public void add(String str, boolean flushIfFull);

  /**
   * appends the character c to the current prefix in O(alphabet size) time. 
   * This method doesn't modify the codebook.
   * @param c: the character to append
   * @return true if the current prefix after appending c is a word in
   * the codebook and false otherwise
   */
    public boolean advance(char c);

    
  /**
   * adds the current prefix to the codebook (if not already a word)
   * then resets the current prefix to the empty string. The running time is 
   * O(1) time.
   * @param flushIfFull a boolean indicating whether the codebook should be 
   *                    flushed when it is at its maximum capacity
   */
    public void add(boolean flushIfFull);

  /**
   * retrieves the codeword corresponding to the current prefix in O(1) time. 
   * Useful in LZW compression.
   * @return the codeword corresponding to the current prefix
   */
  public int getCodeWord();

}