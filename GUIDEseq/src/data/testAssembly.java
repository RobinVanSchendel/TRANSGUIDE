package data;

public class testAssembly {

	public static void main(String[] args) {
		
		String s  =  "";
		System.out.println(s.substring(0, 0));
		System.out.println(s.substring(0, 2));
		
		
		
		Fragment f1 = new Fragment("CTTGAGCTTGGATCAGATTGTCGTTTCCCGCCTTCAGTTTAAACTATCAGTGTTTGACAGTATCTAAAGCAAGGAGAAGCAGAAAGCATAGTGAGGGAATGGTCAGAAGCTCTAGGAGGTGCAGGGATCGATTATCTTGAGAACAAAAGCTGCAATTTTGGAAG");
		Fragment f2 = new Fragment("CTTGAGCTTGGATCAGATTGTCGTTTCCCGCCTTCAGTTTAAACTATCAGTGTTTGACAGTATCTAAAGCAAGGAGAAGCAGAAAGCATAGTGAGGGAATGGTCAGAAGCTCTAGGAGGTGCAGGGATCGATTATCTTGAGAACAAAAGCTGCAATTTTGGAAG");
		Assembler a = new Assembler();
		a.addFragment(f1);
		a.addFragment(f2);
		//System.out.println(f1.length());
		System.out.println(a.toString());
		a.assembleAll();
		System.out.println("Assembled into");
		System.out.println(a.toString());

	}

}
