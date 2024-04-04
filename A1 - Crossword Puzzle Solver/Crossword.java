import java.util.HashMap;
import java.util.Map;

public class Crossword {

 /*
  * -----------------------------------------------------------------------------
  * 
  * Global variables for Crossword class.
  * 
  * -----------------------------------------------------------------------------
  */

 private StringBuilder[] rowStr; // Collection of StringBuilders for K rows
 private StringBuilder[] colStr; // Collection of StringBuilders for K columns
 private int row = 0; // Indexes the current row
 private int col = 0; // Indexes the current column
 private int strValueRow; // Value to determine if row string is a prefix/word/both/neither
 private int strValueCol; // Value to determine if col string is a prefix/word/both/neither
 private boolean isBoardFilledCorrectly; // Flag when board has successfully filled
 private int[] values; // Array of letter values: [0] is the amount the letter appears, [1] is the max amount of times the letter can appear
 private Map<Character, int[]> letterMap = new HashMap<>(); // The map of all characters to determine their frequency
   

 /* Constructuor to declare and initialize structures initially */
 public Crossword() {
  // Sets null for initialization purposes
  rowStr = null;
  colStr = null;
  letterMap = new HashMap<>();
  values = null;
 }

 /* ------------------------------------------------------------------------------------------------------------------
  *                                          fillPuzzle()
  * Checks the current character at the square:
  * - If we have a '+', make sure the placed letter will be a valid placement
  * - If we have a letter, make sure the letter will be a valid placement
  * - If we have a '-', make sure the current board state is valid
  * - If we have a digit, make sure the current letter doesn't appear more than the digit size
  * 
  * returns: 
  *         - The board if we have valid word placements for each row & column
  *         - NULL if there is no valid combinations of words
  * ------------------------------------------------------------------------------------------------------------------*/ 

 public char[][] fillPuzzle(char[][] board, DictInterface dictionary) {
  initializeStructures(board); // Initialize the StringBuilders and map structures once

  if (board[row][col] == '+') {

   /* ------------------------------------------------------------------------------------------------------------------
    *                                          ** OPEN SQUARE '+'**
    * - Adds to the board, stringbuilders, and map
    * - If letter is valid, move to the next square
    * - If letter is invalid, remove it and try the next letter
    * - If every letter was exhausted, return null
    *
    * ------------------------------------------------------------------------------------------------------------------*/

   char letter = 'a'; // Base letter (increments to b,c,... as our search returns unsuccessful)
       
   do {

    // Add current letter to board and StringBuilders
    addLetterToBoard(board, letter); 

    // Get the current letter's frequency of appearance & maximum limit
    values = letterMap.get(letter); 

    if (((values[0] + 1 <= values[1]) || values[1] == 0) && isBelowValid(board, dictionary)) {
     // Appearances did NOT exceed its limit AND square below is valid 

     if (!isLetterValid(board, dictionary)) {
      // Letter is invalid - delete it from the board and StringBuilders
      deleteLetterFromBoard(board, letter);

     } else {
      // Letter is valid - incremenet the letter's appearance in the Map

      if (values[0] < 0) {
       // Value is a garbage value - recompute the frequency of appearance
       values[0] = countLetter(board, letter);
      } else {
       // Value was valid - increment its appearance value
       values[0]++;
      }

      // Save the changes made to the letter's values
      letterMap.put(letter, values);

      if ((!isLastRow(board) && !isLastColumn(board)) || (isLastRow(board) && !isLastColumn(board))) {
       // Go to [row][nextColumn]
       col++;
       fillPuzzle(board, dictionary);
       col--; // Go back a column once we backtrack to this call
      } else if (!isLastRow(board) && isLastColumn(board)) {
       // Go to [nextRow][firstColumn]
       row++;
       int prevCol = col; // store the column temporarily for backtracking purposes
       col = 0;
       fillPuzzle(board, dictionary);
       row--; // Go back a row once we backtrack to this call
       col = prevCol; // Go to the last column once we backtrack to this call

       if (isBoardFilledCorrectly) {
        // Board is filled from backtracking; keep returning the solved puzzle
        return board;
       }
      } else {
       if ((strValueRow == 2 || strValueRow == 3) && (strValueCol == 2 || strValueCol == 3)) {
        // We have valid words for both StringBuilders - return the board
        if (Done(board)) {
         return board;
        }
       }
      }

      if (isBoardFilledCorrectly) {
       // This will continuously return the board each recursive pop off the stack
       return board;
      }

      // We backtracked, and the board isn't finished - delete the current letter
      deleteLetterFromBoard(board, letter);
      values = letterMap.get(letter);
      if (values[0] > 0) {
       // Value was valid - decrement its appearance value
       values[0]--;
      } else {
       // Value is a garbage value - recompute the frequency of appearance
       values[0] = countLetter(board, letter);
      }

      // Save the changes made to the letter's values
      values = letterMap.put(letter, values);
     }
    } else {
     // Square below is invalid: delete it
     deleteLetterFromBoard(board, letter);
    }

   letter++; // Go to next letter

   } while (letter <= 'z');
  } else if (Character.isLetter(board[row][col])) {

   /* ------------------------------------------------------------------------------------------------------------------
    *                                        **PREDEFINED LETTER**
    * - Adds to the StringBuilders, and map
    * - If letter is valid, move to the next square
    * - If letter is invalid, remove it return null
    *
    * ------------------------------------------------------------------------------------------------------------------*/

   // Add the predefined letter to the board
   char predefinedLetter = board[row][col];
   addLetterToBoard(board, predefinedLetter);

   // Get the current letter's frequency of appearance & maximum limit
   values = letterMap.get(board[row][col]);

   if (((values[0] + 1 <= values[1]) || values[1] == 0) && isBelowValid(board, dictionary)) {
    // Appearances is not more than the letter's limit, and square below is valid

    if (!isLetterValid(board, dictionary)) {
      // Letter doesn't work; delete it from the StringBuilders
      rowStr[row].deleteCharAt(col);
      colStr[col].deleteCharAt(row);
    } else {
     if (values[0] < 0) {
      // Value is a garbage value - recompute the frequency of appearance
      values[0] = countLetter(board, predefinedLetter);
     } else {
      // Value was valid - increment its appearance value
      values[0]++;
     }

     // Save the changes made to the letter's values
     letterMap.put(board[row][col], values);

     if ((!isLastRow(board) && !isLastColumn(board)) || (isLastRow(board) && !isLastColumn(board))) {
      // Go to [row][nextColumn]
      col++;
      fillPuzzle(board, dictionary);
      col--; // Go back a column once we backtrack to this call
     } else if (!isLastRow(board) && isLastColumn(board)) {
      // Go to [nextRow][firstColumn]
      row++;
      int prevCol = col;
      col = 0;
      fillPuzzle(board, dictionary);
      row--; // Go back a row once we backtrack to this call
      col = prevCol; // Go to the last column once we backtrack to this call
     } else {
      if ((strValueRow == 2 || strValueRow == 3) && (strValueCol == 2 || strValueCol == 3)) {
       // We have valid words for both StringBuilders - return the board
       if (Done(board)) {
        return board;
       }
      }
     }

     if (isBoardFilledCorrectly) {
      return board;
     }
      // We backtracked, and board isn't done - delete the letter & adjust frequency values
      rowStr[row].deleteCharAt(col);
      colStr[col].deleteCharAt(row);

      // Get the predefined letter's map values
      values = letterMap.get(board[row][col]);

      if (values[0] > 0) {
       // Value was valid - decrement its appearance value
       values[0]--;
      } else {
       // Value is a garbage value - recompute the frequency of appearance
       values[0] = countLetter(board, predefinedLetter);
      }

      // Save the changes made to the letter's values
      values = letterMap.put(board[row][col], values);
    }
   } else {
    // Below isn't valid
    rowStr[row].deleteCharAt(col);
    colStr[col].deleteCharAt(row);
   }
  } else if (board[row][col] == '-') {

   /* ------------------------------------------------------------------------------------------------------------------
    *                                        **UNFILLABLE SQUARE '-'**
    * - Adds to the StringBuilders, and map
    * - If both row and column are valid, go to the next square
    * - If either row or column are invalid, return null
    *
    * ------------------------------------------------------------------------------------------------------------------*/
   
   // Add the '-' to the board and StringBuilders
   char minus = board[row][col];
   addLetterToBoard(board, minus);
   
   // Get the current row's & column's StringBuilders
   StringBuilder rowBuilder = rowStr[row];
   StringBuilder colBuilder = colStr[col];

   if (row == 0) {
    // If '-' is in the first row, automatically make column a word
    strValueCol = 2;
   } else if (!isColWord(board, dictionary, colBuilder)) {
    rowStr[row].deleteCharAt(col);
    colStr[col].deleteCharAt(row);
   }

   if (col == 0) {
    // if '-' is in the first column, automatically make row a word
    strValueRow = 2;
   } else if (!isRowWord(board, dictionary, rowBuilder)) {
    // RowStr[i] + letter is NOT a word -
    if (rowStr[row].length() - 1 >= col) {
     // Delete the '-' from the StringBuilders if we haven't already
     rowStr[row].deleteCharAt(col);
     colStr[col].deleteCharAt(row);
    }
   }

   // Must have row == word, col == word
   if ((strValueRow == 2 || strValueRow == 3) && (strValueCol == 2 || strValueCol == 3)) {
    // We have valid words for both StringBuilders - go to next square or return the
    // board if we are done
    if ((!isLastRow(board) && !isLastColumn(board)) || (isLastRow(board) && !isLastColumn(board))) {
     // We have valid prefixes for both StringBuilders - go to [row][nextColumn]
     col++;
     fillPuzzle(board, dictionary);
     col--; // Go back a column once we backtrack to this call
    } else if (!isLastRow(board) && isLastColumn(board)) {
     // [nextRow][firstColumn]
     row++;
     int prevCol = col;
     col = 0;
     fillPuzzle(board, dictionary);
     row--; // Go back a row once we backtrack to this call
     col = prevCol; // Go to the last column once we backtrack to this call
    } else {
     if (Done(board)) {
      return board;
     }
    }

    if (isBoardFilledCorrectly) {
     // Board is filled from backtracking; keep returning the solved puzzle
     return board;
    }

    // We backtracked, so delete '-'
    rowStr[row].deleteCharAt(col);
    colStr[col].deleteCharAt(row);
   }
  } else if (Character.isDigit(board[row][col])) {


   /* ------------------------------------------------------------------------------------------------------------------
    *                                                **DIGIT**
    * - Adds letter to the board and StringBuilders
    * - If both row and column are valid, and our letter does not exceed the digit's count, go to the next square
    * - If this is false, go to the next letter
    * - If every letter was exhausted, return null
    * ------------------------------------------------------------------------------------------------------------------*/

   // Square is a predefined number
   int digit = board[row][col] - '0';
   char letter = 'a';

   do {
    // Add the letter to the board and StringBuilders
    addLetterToBoard(board, letter);

    // Get the letter's map values, and determine what its limit is
    boolean hasLimit = false;
    values = letterMap.get(letter);
    if (values[1] != 0 && values[0] > 0) {
     // Current letter already has a limit attached to it
     hasLimit = true;
    }

    if ((values[0] + 1 <= digit) && !hasLimit && isBelowValid(board, dictionary)) {
     // Square below is valid
     if (!isLetterValid(board, dictionary)) {
      // Letter is invalid
      deleteLetterFromBoard(board, letter);
     } else {
      // Letter is valid
      if (values[0] < 0) {
       // Value is a garbage value - recompute the frequency of appearance
       values[0] = countLetter(board, letter);
      } else {
       // Value was valid - increment its appearance value
       values[0]++;
      }

      // Digit is letter's new limit
      values[1] = digit;

      // Save the changes
      letterMap.put(letter, values);

      if ((!isLastRow(board) && !isLastColumn(board))) {
       // Go to [row][nextColumn]
       col++;
       fillPuzzle(board, dictionary);
       col--; // Go back a column once we backtrack to this call
      } else if (isLastRow(board) && !isLastColumn(board)) {
       // Go to [row][nextColumn]
       col++;
       fillPuzzle(board, dictionary);
       col--; // Go back a column once we backtrack to this call
      } else if (!isLastRow(board) && isLastColumn(board)) {
       // Go to [nextRow][firstColumn]
       row++;
       int prevCol = col;
       col = 0;
       fillPuzzle(board, dictionary);
       row--; // Go back a row once we backtrack to this call
       col = prevCol; // Go to the last column once we backtrack to this call

       if (isBoardFilledCorrectly) {
        // Board is filled from backtracking; keep returning the solved puzzle
        return board;
       }
      } else {
       values = letterMap.get(letter);
       if ((strValueRow == 2 || strValueRow == 3) && (strValueCol == 2 || strValueCol == 3)) {
        // We have valid words for both StringBuilders - return the board
        if (Done(board)) {
         return board;
        }
       }
      }
     if (isBoardFilledCorrectly) {
      // This will continuously return the board each recursive pop off the stack
      return board;

     }

     // We backtracked, and the board isn't finished - delete the current letter
     deleteLetterFromBoard(board, letter);
     values = letterMap.get(letter);
     if (values[0] > 0) {
      // Value was valid - decrement its appearance value
      values[0]--;
     } else {
      // Value is a garbage value - recompute the frequency of appearance
      values[0] = countLetter(board, letter);
     }
     
     // Reset the letter's limit to 0
     values[1] = 0;
     
     // Save the changes
     values = letterMap.put(letter, values);
     hasLimit = false;
     }
    } else {
     // Square below is invalid: delete it
     deleteLetterFromBoard(board, letter);
    }

    letter++; // Go to next letter

   } while (letter <= 'z');
   board[row][col] = (char) (digit + '0');
  }
   // Nothng valid found
   return null;
 }

  /*
  * ------------------------------------------------------------------------------------------------------------------
  *                                           **checkPuzzle()**
  * Checks each row and column if it's a word
  * - If we have a letter, add it and go to the next square
  * - If we have a '-', stop adding, make sure the current string is a word
  *    - If true, continue going through the row/column
  *    - If false, return false
  * returns:
  * - True if each row and column was a valid word
  * - False if any row or column was not a word
  * ------------------------------------------------------------------------------------------------------------------*/

  public boolean checkPuzzle(char[][] emptyBoard, char[][] filledBoard, DictInterface dictionary) {

   row = 0; // Start at first row
   col = 0; // Start at first column
   StringBuilder checkString = new StringBuilder(""); // Used to hold temporary word
   int strVal; // Value returned from searching the dictionary

   // Add each letter and '-' from each row, check if it's a valid word, if not, return false
   for (int i = 0; i < filledBoard.length; i++) {
    for (int j = 0; j < filledBoard.length; j++) {
     if (emptyBoard[i][j] != '-' && !Character.isDigit(emptyBoard[i][j])) {
      // Block was a fillable block - add it
      checkString.append(filledBoard[i][j]);
     } else if (Character.isDigit(emptyBoard[i][j])) {
      // Block was a digit
      char restrictedLetter = filledBoard[i][j];
      int maximumOccurences = emptyBoard[i][j] - '0';
      if (doesLetterExceedLimit(emptyBoard, filledBoard, restrictedLetter, maximumOccurences)) {
       // Letter exceeds the number at current block
       return false;
      } else {
       // Letter does NOT exceed the number limit at current block, add it to the
       // string
       checkString.append(restrictedLetter);
      }
     } else {
      // Block was unfillable block - make sure we have a word now
      if (checkString.length() != 0) {
       // '-' again; leave
       strVal = dictionary.searchPrefix(checkString, 0, checkString.length() - 1);
       if (strVal != 2 && strVal != 3) {
        // Current string isn't a word - return false
        return false;
       }
      }
      // Current string was a word - empty it and continue
      checkString.setLength(0);
     }
    }

    // Did not run into any '-': check if row is valid word
    if (checkString.length() > 0) {
     strVal = dictionary.searchPrefix(checkString, 0, checkString.length() - 1);
     if (strVal != 2 && strVal != 3) {
      return false;
     }
    }
    // Row was valid
    checkString.setLength(0);
   }

   // Add each letter and '-' from each column, check if it's a valid word, if not, return false
   for (int i = 0; i < filledBoard.length; i++) {
    // Loop through every row index in column
    for (int j = 0; j < filledBoard.length; j++) {
     if (emptyBoard[j][i] != '-') {
      // Block was a fillable block - add it
      checkString.append(filledBoard[j][i]);
     } else {
      // Block was '-' - check if current string is a word
      if (checkString.length() != 0) {
       // '-' doesn't appear consecutively
       strVal = dictionary.searchPrefix(checkString, 0, checkString.length() - 1);
       if (strVal != 2 && strVal != 3) {
        // Current string isn't a word - return false
        return false;
       }
      }

      // Current string was a word - empty it and continue
      checkString.setLength(0);
     }
    }

    // Did not run into any '-' - check if col is valid word
    if (checkString.length() > 0) {
     strVal = dictionary.searchPrefix(checkString, 0, checkString.length() - 1);
     if (strVal != 2 && strVal != 3) {
      return false;
     }
    }
    checkString.setLength(0);
   }
   return true;
  }
  
  /*
   * ------------------------------------------------------------------------------------------------------------------
   * 
   *                          **Helper methods for fillPuzzle() and checkPuzzle()**
   * 
   * ------------------------------------------------------------------------------------------------------------------*/

  /* Initializes the board board and StringBuilders */
  private void initializeStructures(char[][] board) {
   if (rowStr == null || colStr == null) {
    // Runs once at the beginning to initialize structures
    if (rowStr == null) {
     // No array of StringBuilders for rows yet - create the arrays
     rowStr = new StringBuilder[board.length];
     for (int i = 0; i < board.length; i++) {
      rowStr[i] = new StringBuilder();
      for (int j = 0; j < board.length; j++) {
      }
     }
    }

    if (colStr == null) {
     // No arrays of StringBuilders for columns yet - create the arrays
     colStr = new StringBuilder[board.length];
     for (int j = 0; j < board.length; j++) {
      colStr[j] = new StringBuilder();
      for (int i = 0; i < board.length; i++) {
      }
     }
    }
   }

   if (values == null) {
    char[] alphabet = "abcdefghijklmnopqrstuvwxyz".toCharArray();
    for (char characters : alphabet) {
     int[] data = { 0, 0 };
     letterMap.put(characters, data);
    }
   }
  }

  /* Returns true if square below will make a valid word/prefix; false otherwise */
  private boolean isBelowValid(char[][] board, DictInterface dictionary) {
   StringBuilder strBuilder = colStr[col];

   if (isLastRow(board)) {
    // We're on the last row, no need to check anything
    return true;

   } else if (row < board.length - 1 && board[row + 1][col] == '+') {
    // Below is a '+', nothing special needs to be done
    return true;

   } else if (row < board.length - 1 && Character.isLetter(board[row + 1][col])) {
    // Below is a predefined letter - add it to colStr[i]
    colStr[col].append(board[row + 1][col]);

    if (row + 1 < board.length - 1) {
     // Below is NOT the last row - check for valid ColStr[j] prefix or word so far
     if ((row + 2 < board.length - 1) && board[row + 2][col] == '-') {
      // There is a '-' 2 rows down, now we need to check if ColStr[j] is a word
      if (isColWord(board, dictionary, strBuilder)) {
       // The column was a valid word
       colStr[col].deleteCharAt(colStr[col].length() - 1);
       return true;
      }
     } else {
      // There was NOT a '-' 2 rows below, we need to check if ColStr[j] is a prefix
      if (isColPrefix(board, dictionary, strBuilder)) {
       // The column was a valid prefix
       colStr[col].deleteCharAt(colStr[col].length() - 1);
       return true;
      }
     }
    } else {
     // Below is last row - check for valid ColStr[j] word
     if (isColWord(board, dictionary, strBuilder)) {
      // The column was a valid word
      colStr[col].deleteCharAt(colStr[col].length() - 1);
      return true;
     }
    }
   } else if (row < board.length - 1 && board[row + 1][col] == '-') {
    // Block below is a '-' so make sure current column is a valid word
    if (isColWord(board, dictionary, strBuilder)) {
     // The column was a valid word
     return true;
    } else {
     // The column is not a valid word
     return false;
    }
   } else if (row < board.length - 1 && Character.isDigit(board[row + 1][col])) {
    // Square below is a digit, don't do anything
    return true;
   }
   // Added to colStr but wasn't valid; delete the letter added
   colStr[col].deleteCharAt(colStr[col].length() - 1);
   return false;
  }

  /* Adds a letter to the board and StringBuilders */
  private void addLetterToBoard(char[][] board, char letter) {
   // Add to board and StringBuilders
   board[row][col] = letter;
   rowStr[row].append(letter);
   if (colStr[col].length() <= row) {
    // We didn't add the current letter to the colStr yet
    colStr[col].append(letter);
   }
  }

  /* Deletes a letter from the board and StringBuilders */
  private void deleteLetterFromBoard(char[][] board, char letter) {
   board[row][col] = '+';
   rowStr[row].deleteCharAt(col);
   colStr[col].deleteCharAt(row);
  }

  /* Returns true if current square is in last row; false otherwise */
  private boolean isLastRow(char[][] board) {
   if (row < board.length - 1) {
    // Not the last row
    return false;
   }
   // Is the last row
   return true;
  }

  /* Returns true if current square is in last column; false otherwise */
  private boolean isLastColumn(char[][] board) {
   if (col < board[row].length - 1) {
    // Not the last column
    return false;
   }
   // Is the last column
   return true;
  }

  /* Returns true if StringBuilder only has '-'; false othewise */
  private boolean stringIsAllUnfillableSquares(StringBuilder strBuilder) {
   for (int i = 0; i < strBuilder.length(); i++) {
    // Search the entire StringBuilder
    if (strBuilder.charAt(i) != '-') {
     // The character is not a '-'
     return false;
    }
   }
   // StringBuilder is filled completely with '-'
   return true;
  }

  /* Returns true if colStr[j] + letter is a word; false otherwise */
  private boolean isColWord(char[][] board, DictInterface dictionary, StringBuilder colBuilder) {
   int lastUnfillableIndex = findLastUnfillableBlock(colBuilder); // Check the last '-' index

   if (stringIsAllUnfillableSquares(colBuilder) || (((row > 0) && (board[row - 1][col] == '-'))
     && ((row < board.length - 1) && (board[row + 1][col] == '+')))) {
    // Current column is only '-', OR there exists a '-' right before the current square and the square after is fillable only
    strValueCol = 3;
    return true; // ColStr is automatically a word and prefix

   } else if (board[row][col] != '-' && lastUnfillableIndex == -1) {
    // Current square is a letter, and StringBuilder has no '-' before it
    strValueCol = dictionary.searchPrefix(colBuilder, 0, colBuilder.length() - 1); // [start -> end]

   } else if ((board[row][col] != '-') && (lastUnfillableIndex != -1)) {
    // Current square is a letter, and a '-' exists somewhere before current square
    // in the StringBuilder
    strValueCol = dictionary.searchPrefix(colBuilder, lastUnfillableIndex + 1, colBuilder.length() - 1); // (last '-' -> end]

   } else if ((board[row][col] == '-') && (lastUnfillableIndex == -1)) {
    // Current square is a '-', and StringBuilder has no '-' before it
    strValueCol = dictionary.searchPrefix(colBuilder, 0, colBuilder.length() - 2); // [start -> end)

   } else if (board[row][col] == '-' && (row > 0) && board[row - 1][col] == '-') {
    // Current square is a '-', and '-' right before it, so colStr was a word
    // already
    strValueCol = 3;
    return true; // ColStr is automatically a word and prefix

   } else if ((board[row][col] == '-') && (lastUnfillableIndex != -1)) {
    // Current square is a '-', and a '-' exists somewhere before current square in
    // the StringBuilder
    strValueCol = dictionary.searchPrefix(colBuilder, lastUnfillableIndex + 1, colBuilder.length() - 2); // (last '-' -> end)
   }

   if (strValueCol != 2 && strValueCol != 3) {
    // Current colStr is NOT a word
    return false;
   }
   // Current colStr is a word
   return true;
  }

  /* Returns true if colStr[j] + letter is a prefix; false otherwise */
  private boolean isColPrefix(char[][] board, DictInterface dictionary, StringBuilder colBuilder) {
   int lastUnfillableIndex = findLastUnfillableBlock(colBuilder); // Check the last '-' index

   if (stringIsAllUnfillableSquares(colBuilder) || (((row > 0) && (board[row - 1][col] == '-'))
     && ((row < board.length - 1) && (board[row + 1][col] == '+')))) {
    // Current column is only '-', OR their exists a '-' right before the current square and not after
    strValueCol = 3;
    return true; // colStr is automatically a word and prefix

   } else if (board[row][col] != '-' && lastUnfillableIndex == -1) {
    // Current block has a letter, and no unfillable blocks in the column
    strValueCol = dictionary.searchPrefix(colBuilder, 0, colBuilder.length() - 1); // [start -> end]

   } else {
    // Current block is a letter and has an unfillable block in the column
    strValueCol = dictionary.searchPrefix(colBuilder, lastUnfillableIndex + 1, colBuilder.length() - 1); // (last '-' -> end]
   }

   if (strValueCol != 3 && strValueCol != 1) {
    // strCol is NOT a prefix
    return false;
   }
   // strCol is a prefix
   return true;
  }

  /* Returns true if rowStr[i] + letter is a word; false otherwise */
  private boolean isRowWord(char[][] board, DictInterface dictionary, StringBuilder rowBuilder) {
   int lastUnfillableIndex = findLastUnfillableBlock(rowBuilder); // Check the last '-' index

   if (stringIsAllUnfillableSquares(rowBuilder) || (((col > 0) && (board[row][col - 1] == '-'))
     && ((col < board.length - 1) && (board[row][col + 1] != '-')))) {
    // Current row is only '-', or their exists a '-' right before the currentsquare and not after
    strValueRow = 3;
    return true; // rowStr is automatically a word and prefix

   } else if (board[row][col] != '-' && lastUnfillableIndex == -1) {
    // Current square is a letter, and StringBuilder has no '-' before it
    strValueRow = dictionary.searchPrefix(rowBuilder, 0, rowBuilder.length() - 1); // [start -> end]

   } else if ((board[row][col] != '-') && (lastUnfillableIndex != -1)) {
    // Current square is a letter, and a '-' exists somewhere before current square
    // in the StringBuilder
    strValueRow = dictionary.searchPrefix(rowBuilder, lastUnfillableIndex + 1, rowBuilder.length() - 1); // (last '-' -> end]

   } else if ((board[row][col] == '-') && (lastUnfillableIndex == -1)) {
    // Current square is a '-', and StringBuilder has no '-' before it
    strValueRow = dictionary.searchPrefix(rowBuilder, 0, rowBuilder.length() - 2); // [0 -> end)

   } else if ((board[row][col] == '-') && ((col > 0) && (board[row][col - 1] == '-'))) {
    // Current square is a '-', and '-' right before it, so rowStr was a word
    // already
    strValueRow = 3;
    return true; // rowStr is automatically a word and prefix

   } else if ((board[row][col] == '-') && (lastUnfillableIndex != -1)) {
    // Current square is a '-', and a '-' exists somewhere before current square in
    // the StringBuilder
    strValueRow = dictionary.searchPrefix(rowBuilder, lastUnfillableIndex + 1, rowBuilder.length() - 2); // (last -> end)

   }
   if (strValueRow != 2 && strValueRow != 3) {
    // rowStr is NOT a word
    return false;
   }
   // rowStr is a word
   return true;
  }

  /* Returns true if rowStr[i] + letter is a prefix; false otherwise */
  private boolean isRowPrefix(char[][] board, DictInterface dictionary, StringBuilder rowBuilder) {
   int lastUnfillableIndex = findLastUnfillableBlock(rowBuilder); // Check the last '-' index

   if (stringIsAllUnfillableSquares(rowBuilder) || (((col > 0) && (board[row][col - 1] == '-'))
     && ((col < board.length - 1) && (board[row][col + 1] != '-')))) {
    // Current row is only '-', or their exists a '-' right before the currentsquare
    // and not after
    strValueRow = 3;
    return true; // rowStr is automatically a word and prefix

   } else if (board[row][col] != '-' && lastUnfillableIndex == -1) {
    // Current square is a letter, and no '-' in the StringBuilder
    strValueRow = dictionary.searchPrefix(rowBuilder, 0, rowBuilder.length() - 1); // [start -> end]

   } else {
    // Current block is a letter and has an unfillable block in the column
    strValueRow = dictionary.searchPrefix(rowBuilder, lastUnfillableIndex + 1, rowBuilder.length() - 1); // (last '-' -> end]
   }

   if (strValueRow != 3 && strValueRow != 1) {
    // rowStr is NOT a prefix
    return false;
   }
   // rowStr is a prefix
   return true;
  }

  /* Returns true if there are valid words in both rows and column; false otherwise */
  private boolean Done(char[][] board) {
   if (row == board.length - 1 && col == board.length - 1) {
    // We've reached the solution; return the board
    isBoardFilledCorrectly = true;
    return true;
   }
   // We've yet to reach the end
   return false;
  }

  /* Returns true if a letter exceeds the maximum amount of times it's allowed to appear; false otherwise */
  private boolean doesLetterExceedLimit(char[][] emptyBoard, char[][] filledBoard, char letter, int maximumOccurences) {
   int currentOccurences = 0; // Keeps track of how many times the restricted letter occurs on the board

   // Search the board
   for (int i = 0; i < filledBoard.length; i++) {
    for (int j = 0; j < filledBoard.length; j++) {
     if (filledBoard[i][j] == letter) {
      currentOccurences++; // Found a match
      if (currentOccurences > maximumOccurences) {
       // Appears too many times in the board - invalid board
       return true;
      }
     }
    }
   }
   // Appears less than or equal to the maximum limit - valid board
   return false;
  }

  /* Checks if the placed letter is valid in both row/col StringBuilders */
  private boolean isLetterValid(char[][] board, DictInterface dictionary) {
   StringBuilder colBuilder = colStr[col]; // StringBuilder of current column
   StringBuilder rowBuilder = rowStr[row]; // StringBuilder of current row

   // Check the column
   if (isLastRow(board)) {
    // Current row is the last row
    if (!isColWord(board, dictionary, colBuilder)) {
     // Current column StringBuilder is NOT a word
     return false;
    }
   } else if ((row < board.length) && board[row + 1][col] == '-') {
    // The current square is NOT the last row but there is a '-' below it
    if (!isColWord(board, dictionary, colBuilder)) {
     // Current column Stringbuilder is NOT a word
     return false;
    }
   } else {
    // Current row is NOT the last row
    if (!isColPrefix(board, dictionary, colBuilder)) {
     // Current column StringBuilder is NOT a prefix
     return false;
    }
   }
   
   // Check the row
   if (isLastColumn(board)) {
    // Current column is the last column
    if (!isRowWord(board, dictionary, rowBuilder)) {
     // Current row StringBuidler is NOT a word
     return false;
    }
   } else if ((col < board.length) && board[row][col + 1] == '-') {
    // The current square is NOT in the last column but there is a '-' to the right
    // of it
    if (!isRowWord(board, dictionary, rowBuilder)) {
     // Current row StringBuilder is NOT a prefix
     return false;
    }
   } else {
    // Current column is NOT the last column
    if (!isRowPrefix(board, dictionary, rowBuilder)) {
     // Current row StringBuilder is NOT a prefix
     return false;
    }
   }
   // Letter is valid
   return true;
  }

  /* If there are any '-' blocks before current square; finds the closest and returns its position within the StringBuilder */
  private int findLastUnfillableBlock(StringBuilder strBuilder) {
   int lastUnfillableIndex;
   for (int i = strBuilder.length() - 2; i >= 0; i--) {
    if ((strBuilder.charAt(i) == '-')) {
     // '-' is in the StringBuilder and not the current square
     lastUnfillableIndex = i; // Update where the latest '-' is
     return lastUnfillableIndex;
    }
   }
   // Did not find any '-'
   return -1;
  }

  /* Counts the amount of times a letter occurs on the board */
  private int countLetter(char[][] board, char letter) {
   int count = 0;
   for (int i = 0; i < board.length; i++) {
    for (int j = 0; j < board.length; j++) {
     if (letter == board[i][j]) {
      count++; // Found a match
     }
    }
   }
   return count;
  }
}