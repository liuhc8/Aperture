
public class CombinationCalculator {
	
	private static long factorial(int start,int end) {
		long sum=1L;
		for(int i=start;i<=end;++i) {
			sum*=i;
		}
		return sum;
	}
	
	public static long combination(int m,int n) {
		return m <= n ? factorial(n-m+1,n) / factorial(1,m) : 0;
	}
	
}
