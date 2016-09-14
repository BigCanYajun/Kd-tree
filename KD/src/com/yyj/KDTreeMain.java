package com.yyj;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class KDTreeMain {
	public static int KDTCount = 0; // ͳ����kdt ������ʱ�򣬼����˺ͼ�����ľ���
	public static void main(String[] args) {
		/*
		 *  n >> 2^xn ʱ�� KDTCount������ < n
		 *  
		 *	=========================
		 *	n = 40000, xn = 10
			buld kdt time = 30760.0
			
			query kdt time = 1404.0
			best  = 0.3296984744447501
			KDTCount = 4488
			
			query brute time = 3317.0
			best2 = 0.3296984744447501
	
			==========================
			n = 50000, xn = 10
			buld kdt time = 49664.0
			
			query kdt time = 558.0
			best  = 0.3355435846472523
			KDTCount = 2056
			
			query brute time = 5557.0
			best2 = 0.3355435846472523
			
			==========================
			n = 50000. xn = 20
			buld kdt time = 63560.0
			
			query kdt time = 15136.0
			best  = 0.8319764077450744
			KDTCount = 37500
			
			query brute time = 5791.0
			best2 = 0.8319764077450744

		 **/
		int n = 50000;  // ���������
		int xn = 10;	// ������ά��
		int deep = 0;	// ��
	
		// �������ѵ����������
		List<Point> pointList = new LinkedList<Point>();
		for (int i = 0; i < n; i++) {
			double[] d = new double[xn];
			for (int j = 0; j < d.length; j++) {
				d[j] = Math.random();
			}
			pointList.add(new Point(d));
		}
		
		// build tree
		System.out.println("beging insert...");
		double t1 = System.currentTimeMillis();
		KDTreeMain kdt = new KDTreeMain();
		Node root = new Node();
		kdt.insert(root, pointList, deep);
		double t2 = System.currentTimeMillis();
		System.out.println("buld kdt time = " + (t2 - t1));
		
		// show tree
//		char[] path = new char[30];
//		int pi = 0;
//		showKDTree(root, path, pi);
		
		// Ŀ���
		double[] f = new double[xn];
		for (int j = 0; j < f.length; j++) {
			f[j] = Math.random();
		}
		Point p = new Point(f);

		// KDT����
		double t3 = System.currentTimeMillis();
		double best = Double.MAX_VALUE;
		best = query(root, p, best, deep);
		double t4 = System.currentTimeMillis();
		System.out.println("\nquery kdt time = " + (t4 - t3));
		System.out.println("best  = " + best);
		System.out.println("KDTCount = " + KDTCount);
		
		// ������
		double t5 = System.currentTimeMillis();
		int index = 0;
		double best2 = Double.MAX_VALUE;
		for (int i = 0; i < n; i++) {
			double dist = getDist(p, pointList.get(i));
			if (dist < best2) {
				best2 = dist;
				index = i;
			}
		}
		double t6 = System.currentTimeMillis();
		System.out.println("\nquery brute time = " + (t6 - t5));
		System.out.println("best2 = " + best2);
		// System.out.println("goal point = " + p.x[0] + " , " + p.x[1]);
		// System.out.println("neast point = " + pointList.get(index).x[0] + " , " + pointList.get(index).x[1]);
	}
	
	// build kdtree
	private void insert(Node root, List<Point> pointList, int deep) {
		int mid = pointList.size() / 2;
		
		// ������õ���λ��
		Point.deep = deep;
		Collections.sort(pointList);
		
		// ���ƿ��ŵķ����õ���λ��
		// getMedian(pointList, 0, pointList.size() - 1, mid, deep);
		// showList(pointList);
		// System.out.println("=========================");
		int pl = mid;
		int pr = mid;
		while(pl >= 0 && pointList.get(pl).x[deep] == pointList.get(mid).x[deep]) pl--;
		while(pr < pointList.size() && pointList.get(pr).x[deep] == pointList.get(mid).x[deep]) pr++;
		List<Point> pointListLeft = pointList.subList(0, pl + 1);
		List<Point> pointListMid = pointList.subList(pl + 1, pr);
		List<Point> pointListRight = pointList.subList(pr, pointList.size());
		
		root.pointList = pointListMid;
		if (pointListLeft.size() > 0) {
			root.l = new Node();
			insert(root.l, pointListLeft, (deep + 1) % pointList.get(0).x.length);
		}
		if (pointListRight.size() > 0) {
			root.r = new Node();
			insert(root.r, pointListRight, (deep + 1) % pointList.get(0).x.length);
		}
		
	}
	
	// search the nearest point to p in KDTree
	private static double query(Node root, Point p, double best, int deep) {
		if (root == null) return Double.MAX_VALUE;  
	    double dist;  
	    if (root.l == null && root.r == null) {  
	        for (int i = 0; i < root.pointList.size(); i++) {  
	            KDTCount++;  
	            dist = getDist(root.pointList.get(i), p);  
	            best = dist < best ? dist : best;  
	        }  
	        return best;  
	    }  
	  
	    // left or right  
	    if (p.x[deep] <= root.pointList.get(0).x[deep]) {  
	        best = query(root.l, p, best, (deep + 1) % p.x.length);
	    } else {  
	        best = query(root.r, p, best, (deep + 1) % p.x.length);
	    }  
	    // cur  
	    for (int i = 0; i < root.pointList.size(); i++) {  
	        KDTCount++;  
	        dist = getDist(root.pointList.get(i), p);  
	        best = dist < best ? dist : best;  
	    }  
	    // another side  
	    if (best >= Math.abs(p.x[deep] - root.pointList.get(0).x[deep])) {  
	        double distAnother = Double.MAX_VALUE;  
	        if (p.x[deep] <= root.pointList.get(0).x[deep]) {  
	            distAnother = query(root.r, p, best, (deep + 1) % p.x.length);
	        } else {  
	            distAnother = query(root.l, p, best, (deep + 1) % p.x.length);
	        }  
	        if (distAnother < best) {  
	            best = distAnother;  
	        }  
	    }  
	    return best;  
	}
	
	// print kdtree
	private static void showKDTree(Node root, char[] path, int pi) {
		if (root == null) return;
		System.out.print(pi + "# ");
		for (int i = 0; i < pi; i++) {
			System.out.print(path[i] + " ");
		}
		// mid
		showList(root.pointList);
		// left
		path[pi++] = 'L';
		showKDTree(root.l, path, pi);
		pi--;
		// right
		path[pi++] = 'R';
		showKDTree(root.r, path, pi);
		pi--;
	}
	// ŷʽ����
	private static double getDist(Point p1, Point p2) {
		double sum = 0;
		for (int i = 0; i < p1.x.length; i++) {
			sum += (p1.x[i] - p2.x[i]) * (p1.x[i] - p2.x[i]);
		}
		if (sum == 0) return Double.MAX_VALUE;
		return Math.sqrt(sum);
	}
	// ���ƿ��ŵ�˼���õ���λ����O(n)ʱ�临�Ӷ�
	private void getMedian(List<Point> pointList, int l, int r, int k, int deep) {
		if (l == r && k == 0) return;  
	    int pl = l;  
	    int pr = r;  
	    double[] tmp = pointList.get(l).x;  
	    while (pl < pr) {  
	        while (pl < pr && pointList.get(pr).x[deep] > tmp[deep]) pr--;  
	        if (pl >= pr) break;  
	        pointList.get(pl++).x = pointList.get(pr).x;  
	        while (pl < pr && pointList.get(pl).x[deep] < tmp[deep]) pl++;  
	        if (pl >= pr) break;  
	        pointList.get(pr--).x = pointList.get(pl).x;
	    }  
	    pointList.get(pl).x = tmp;  
	  
	    if(pl - l == k) return;  
	    if(pl - l >  k) {  
	        getMedian(pointList, l, pl - 1, k, deep);  
	    } else {  
	        getMedian(pointList, pl + 1, r, k - (pl - l + 1), deep);  
	    }  
	}
	// ��ӡһ�����б�
	private static void showList(List<Point> pointList) {
		for (int i = 0; i < pointList.size(); i++) {
			for( int j = 0; j < pointList.get(i).x.length; j++) {
				System.out.print(pointList.get(i).x[j] + ",");
			}
			System.out.print(" / ");
		}
		System.out.println();
	}
}
// kdtree��Ľڵ�
class Node {
	List<Point> pointList = new LinkedList<Point>();
	Node l = null;
	Node r = null;
}
// ���ݵ�
class Point implements Comparable<Point>{
	public static int deep = 0;
	double[] x;
	public Point(double[] d) {
		x = new double[d.length];
		for (int i = 0; i < d.length; i++) {
			x[i] = d[i];
		}
	}
	public int compareTo(Point o) {
		// return (int)(this.x[deep] == other.x[deep]); ������Ϊx��ֵ��0~1֮�䣬��ôint����0��
		Point other = (Point)o;
		if (this.x[deep] == other.x[deep]) return 0;
		if (this.x[deep] >  other.x[deep]) return 1;
		return -1;
	}
}

