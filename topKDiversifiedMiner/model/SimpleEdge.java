package model;

public class SimpleEdge implements Comparable<SimpleEdge> {
	 
	public Integer vlabel1 = 0;
	public Integer vlabel2 = 0;
	@Override
	public int compareTo(SimpleEdge o) {
		// TODO Auto-generated method stub
		int flag = vlabel1.compareTo(o.vlabel1);
		if(flag == 0) {
			return vlabel2 - o.vlabel2;
		} else {
			return flag;
		}
	} 

}
