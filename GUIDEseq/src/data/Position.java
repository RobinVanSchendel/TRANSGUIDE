package data;

public class Position {
	private String chr;
	private int position;
	public Position(String c, int pos) {
		this.chr = c;
		this.position = pos;
	}
	public String getChr() {
		return chr;
	}
	public int getPosition() {
		return position;
	}
}
