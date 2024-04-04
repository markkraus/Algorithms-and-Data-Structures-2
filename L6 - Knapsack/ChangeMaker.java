import java.io.File;
import java.util.Scanner;
import java.util.ArrayList;
public final class ChangeMaker {

	private int[] coins;
	private int nCoins;
	private int change;
	private int[] CM; //the dynamic programming table
	private int[] decision;//remember the decision that we make at each subproblem
	private ArrayList<Integer> selectedCoins;

	public static void main(String[] args){
		if(args.length != 1){
			System.out.println("Usage: java ChangeMaker <filename>");
		} else {
			new ChangeMaker(args[0]);
		}
	}

	public ChangeMaker(String filename){
		Scanner scan = null;
		try {
				scan = new Scanner(new File(filename));
				nCoins = Integer.parseInt(scan.nextLine());
				coins = new int[nCoins];
				for(int i=0; i<nCoins; i++){
					coins[i] = Integer.parseInt(scan.nextLine());
				}
				change = Integer.parseInt(scan.nextLine());
				CM = new int[change+1];
				decision = new int[change+1];
				selectedCoins = new ArrayList<>();

		} catch (Exception e) {
			System.out.println("Error reading file " + filename + ": " + e.getMessage());
		} finally{
      scan.close();
    }

		solve();

		System.out.println("The minimum number of coins to make a change of value "
											+ change + " is " + CM[change]);
		buildSolution();
		System.out.println("The selected coins are: " + selectedCoins.toString());

	}

	private void solve(){
		CM[0] = 0; //0 coins are needed for a change value of zero
		for(int i=1; i<= change; i++){
			int min = Integer.MAX_VALUE;
			for(int j=0; j<nCoins; j++){
				if(i >= coins[j]){ //i is the current change value
					if(1 + CM[i - coins[j]] < min){ //coin j taken
						//TODO: Update min and CM[i]
            min = 1 + CM[i - coins[j]];
            CM[i] = min;
						decision[i] = j;
					}
				}
			}
		}
	}

	private void buildSolution(){
		int cur = change; //start from the largest subproblem
		for(int i=0; i<CM[change]; i++){ //CM[change] is the minimum number of coins
				//decision[cur] is the coin number selected at subproblem cur
				//TODO: Add the selected coin value (coins[decision[cur]]) to selectedCoins
        selectedCoins.add(coins[decision[cur]]);
				cur = cur - coins[decision[cur]]; //the next subproblem is cur - coins[decision[cur]]
		}
	}

}
