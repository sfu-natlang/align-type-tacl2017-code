import java.io.Serializable;
/**
 * This class implements a Triple object
 * @author anahita
 *
 * @param <L>
 * @param <M>
 * @param <R>
 */

public class Triple<L,M,R> implements Serializable{
	final L left;
	final M middle;
	final R right;
	
	public Triple(L left, M middle, R right) {
		super();
		this.left = left;
		this.middle = middle;
		this.right = right;
	}
	/*public void setTriple(L left, M middle, R right){
		this.left = left;
		this.middle = middle;
		this.right = right;
	}*/
	@Override
	public int hashCode(){
		int prime = 31;
		int result = 1;
		result = prime * result + ((left == null) ? 0 : left.hashCode());
		result = prime * result + ((middle == null) ? 0 : middle.hashCode());
		result = prime * result + ((right == null) ? 0 : right.hashCode());
		return result;
	}
	
	@SuppressWarnings("rawtypes")
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Triple other = (Triple) obj;
		
		if (left == null) {
			if (other.left != null)
				return false;
		} else if (!left.equals(other.left))
			return false;
		if (middle == null) {
			if (other.middle != null)
				return false;
		} else if (!middle.equals(other.middle))
			return false;
		if (right == null) {
			if (other.right != null)
				return false;
		} else if (!right.equals(other.right))
			return false;
		return true;
	}
	public String toString(){
		return "(" + left + "," + middle +"," + right +")";
	}
}
