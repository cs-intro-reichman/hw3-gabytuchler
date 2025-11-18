/** Functions for checking if a given string is an anagram. */
public class Anagram {
	public static void main(String args[]) {
		// Tests the isAnagram function.
		System.out.println(isAnagram("silent","listen"));  // true
		System.out.println(isAnagram("William Shakespeare","I am a weakish speller")); // true
		System.out.println(isAnagram("Madam Curie","Radium came")); // true
		System.out.println(isAnagram("Tom Marvolo Riddle","I am Lord Voldemort")); // true

		// Tests the preProcess function.
		System.out.println(preProcess("What? No way!!!"));
		
		// Tests the randomAnagram function.
		System.out.println("silent and " + randomAnagram("silent") + " are anagrams.");
		
		// Performs a stress test of randomAnagram 
		String str = "1234567";
		Boolean pass = true;
		//// 10 can be changed to much larger values, like 1000
		for (int i = 0; i < 10; i++) {
			String randomAnagram = randomAnagram(str);
			System.out.println(randomAnagram);
			pass = pass && isAnagram(str, randomAnagram);
			if (!pass) break;
		}
		System.out.println(pass ? "test passed" : "test Failed");
	}  

	// Returns true if the two given strings are anagrams, false otherwise.
	public static boolean isAnagram(String str1, String str2) {
			str1 = preProcess(str1);
			str2 = preProcess(str2);
			if (str1.length() != str2.length()){
				return false; }
			for (int i=0; i<str1.length();i++){
				int countstr1 = 0;
				int countstr2 = 0;
				 for (int j=0; j<str2.length();j++){
					if (str1.charAt(i) == str2.charAt(j)){
						 countstr2 = countstr2 + 1;}}
				 for (int k=0; k<str1.length();k++){
					if (str1.charAt(i) == str1.charAt(k)){
						 countstr1 = countstr1+ 1;}}
				if (countstr1 == 0 || countstr2 == 0){
					return false;}
				if (countstr1!=countstr2){
					return false;}}
			return true;
		// Replace the following statement with your code
	}
	   
	// Returns a preprocessed version of the given string: all the letter characters are converted
	// to lower-case, and all the other characters are deleted, except for spaces, which are left
	// as is. For example, the string "What? No way!" becomes "whatnoway"
	public static String preProcess(String str) {
		String result = "";
		for (int i=0; i<str.length();i++){
			char c = str.charAt(i);
			if (Character.isLetter(c)){
				if (Character.isUpperCase(c)){
					c = Character.toLowerCase(c);}
				result = result + c;}
			}
		return result;
	} 
	   
	// Returns a random anagram of the given string. The random anagram consists of the same
	// characters as the given string, re-arranged in a random order. 
	public static String randomAnagram(String str) {
		String newstr = "";
		while (str.length()>0){
			int random = (int)(Math.random() * str.length());
			String newchar = ("" + str.charAt(random));
			newstr = newstr + newchar;
			str = str.substring(0, random) + str.substring(random + 1);}
		return newstr;
	}
}
