// Computes the periodical payment necessary to pay a given loan.
public class LoanCalc {
	
	static double epsilon = 0.001;  // Approximation accuracy
	static int iterationCounter;    // Number of iterations 
	
	// Gets the loan data and computes the periodical payment.
    // Expects to get three command-line arguments: loan amount (double),
    // interest rate (double, as a percentage), and number of payments (int).  
	public static void main(String[] args) {		
		// Gets the loan data
		double loan = Double.parseDouble(args[0]);
		double rate = Double.parseDouble(args[1]);
		int n = Integer.parseInt(args[2]);
		System.out.println("Loan = " + loan + ", interest rate = " + rate + "%, periods = " + n);

		// Computes the periodical payment using brute force search
		System.out.print("\nPeriodical payment, using brute force: ");
		System.out.println((int) bruteForceSolver(loan, rate, n, epsilon));
		System.out.println("number of iterations: " + iterationCounter);

		// Computes the periodical payment using bisection search
		System.out.print("\nPeriodical payment, using bi-section search: ");
		System.out.println((int) bisectionSolver(loan, rate, n, epsilon));
		System.out.println("number of iterations: " + iterationCounter);
	}

	// Computes the ending balance of a loan, given the loan amount, the periodical
	// interest rate (as a percentage), the number of periods (n), and the periodical payment.
	private static double endBalance(double loan, double rate, int n, double payment) {	
		
		for (int i=0; i<n; i++){
			loan = (loan - payment) * (1 + rate) ;	
		}
		return loan;
	}
	
	// Uses sequential search to compute an approximation of the periodical payment
	// that will bring the ending balance of a loan close to 0.
	// Given: the sum of the loan, the periodical interest rate (as a percentage),
	// the number of periods (n), and epsilon, the approximation's accuracy
	// Side effect: modifies the class variable iterationCounter.
    public static double bruteForceSolver(double loan, double rate, int n, double epsilon) {
			double g = loan / n;
			LoanCalc.iterationCounter = 0;
			while (endBalance(loan, rate, n, g)>0) {
				g = g + epsilon;
				LoanCalc.iterationCounter = LoanCalc.iterationCounter + 1;
			}
		return g;
    }
    
    // Uses bisection search to compute an approximation of the periodical payment 
	// that will bring the ending balance of a loan close to 0.
	// Given: the sum of the loan, the periodical interest rate (as a percentage),
	// the number of periods (n), and epsilon, the approximation's accuracy
	// Side effect: modifies the class variable iterationCounter.
    public static double bisectionSolver(double loan, double rate, int n, double epsilon) {  
		//double low=loan/n;	
		//double high=(loan*(1+rate))/n;
		//int iterationCounter = 0;
		//double g = (low + high)/ 2;
		//while ((high-low) > epsilon ){
			
		//}
		//return g;
		iterationCounter = 0; 

    // 2. קביעת גבולות הטווח [L, H]
    
    // L (low): גבול תחתון. f(L) > 0 (תשלום הקרן בלבד).
    double low = loan / n; 
    
    // H (high): גבול עליון. f(H) < 0 (תשלום שמבטיח יתרה שלילית).
    // נבחר loan * (1 + rate), כפי שהוסבר, שמבטיח כי f(H) < 0.
    double high = loan * (1.0 + rate);
    
    // 3. משתנה עזר לנקודת האמצע (g)
    double g = 0.0;

    // 4. לולאת החיפוש: ממשיכים כל עוד רוחב הטווח גדול מהדיוק (epsilon).
    while ((high - low) > epsilon) {
        
        // א. עדכון המונה
        iterationCounter++;
        
        // ב. חישוב נקודת האמצע (g)
        g = (low + high) / 2.0; 
        
        // ג. חישוב היתרה בנקודת האמצע (f(g))
        double balanceG = endBalance(loan, rate, n, g);
        
        // 5. עדכון הטווח (חציית הטווח)
        // מכיוון ש-f היא מונוטונית יורדת:
        
        if (balanceG > 0) {
            // אם היתרה חיובית, התשלום נמוך מדי.
            // הפתרון הוא מימין ל-g. מעבירים את L ל-g.
            low = g;
        } else {
            // אם היתרה שלילית או 0, התשלום מספיק.
            // הפתרון הוא משמאל ל-g. מעבירים את H ל-g.
            high = g;
        }
    }
    
    // 6. החזרת הקירוב: הערך g האחרון או הממוצע של הגבולות הסופיים.
    // (הערה: L ו-H כבר קרובים מאוד זה לזה, כך ש-g מייצג את הקירוב הטוב ביותר).
    return g;
    }
}