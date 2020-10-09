import java.util.Arrays;

public class LoserTree{
    private int[] tree,leaves;
    private BpIndexList.GetAndSetItr[] itrList;
    private int size;

    public LoserTree(BpIndexList[] indexMatrix) {
    	int size = indexMatrix.length;
    	this.size=size;
        this.leaves=new int[size];
        this.tree = new int[size];
        this.itrList=new BpIndexList.GetAndSetItr[size];
        
        Arrays.fill(tree, -1);
        
        for (int i = 0; i < size; i++) {
        	itrList[i]=indexMatrix[i].getGetAndSetItr();
        	itrList[i].moveToNext();
        	leaves[i] = itrList[i].getCurrentCode();
        }
        for (int i = size - 1; i >= 0; i--) {
        	adjust(i);
        }
    }


    private void adjust(int s) {
        int t = (s + size) / 2;
        while (t > 0) {
            if (s >= 0 && (tree[t] == -1 || leaves[s]>leaves[tree[t]])) {
                int temp = s;
                s = tree[t];
                tree[t] = temp;
            }
            t /= 2;
        }
        tree[0] = s;
    }


    private void add(int newCode, int s) {
        leaves[s]=newCode;
        adjust(s);
    }


    private void del(int s) {
    	this.size--;
    	int[] newLeaves=new int[size];
    	System.arraycopy(leaves, 0, newLeaves, 0, s);
    	System.arraycopy(leaves, s+1, newLeaves, s, size-s);
    	this.leaves=newLeaves;
    	
    	BpIndexList.GetAndSetItr[] newItrList=new BpIndexList.GetAndSetItr[size];
    	System.arraycopy(itrList, 0, newItrList, 0, s);
    	System.arraycopy(itrList, s+1, newItrList, s, size-s);
    	this.itrList=newItrList;

        this.tree = new int[size];

        Arrays.fill(tree, -1);
        for (int i = size - 1; i >= 0; i--) {
            adjust(i);
        }

    }


    private int getLeafCode(int s) {
        return leaves[s];
    }
    
    private BpIndexList.GetAndSetItr getLeafIter(int s) {
        return itrList[s];
    }


    private int getWinner() {
        return tree.length > 0 ? tree[0] : -1;
    }
    
    public BpIndexList mergeAndSortIndexs() {
    	BpIndexList finalIndexList=new BpIndexList((int)CombinationCalculator.combination(5, 32));
    	
    	long mergedStart=0L;
    	int winner=getWinner();
    	BpIndexList.GetAndSetItr winnerItr=getLeafIter(winner);
    	winnerItr.setCurrentMergedStart(mergedStart);
    	int winnerCode=getLeafCode(winner);
    	//System.out.println("######"+Integer.toHexString(winnerCode)+":::"+winner+":::"+winnerItr.getCurrentInitStart()+":::"+winnerItr.getCurrentLen()+":::"+mergedStart);
    	finalIndexList.addCodeAndStartPos(winnerCode, mergedStart);
    	//System.out.println(Integer.toHexString(winnerCode)+":::"+mergedStart);
    	mergedStart+=winnerItr.getCurrentLen();
    	
    	
    	while (true) {
    		while(true) {
    			if (winnerItr.hasNext()) {
    				winnerItr.moveToNext();
                	int nextCode=winnerItr.getCurrentCode();
                	if(nextCode==winnerCode) {
                		winnerItr.setCurrentMergedStart(mergedStart);
                		//System.out.println("######"+Integer.toHexString(winnerCode)+":::"+winner+":::"+winnerItr.getCurrentInitStart()+":::"+winnerItr.getCurrentLen()+":::"+mergedStart);
                    	mergedStart+=winnerItr.getCurrentLen();
                	}else {
                		add(nextCode,winner);
                		break;
                	}
                } else {
                	del(winner);
                	break;
                }
    		}

            winner = getWinner();
            if (winner==-1) {
            	finalIndexList.setLastLen(mergedStart);
                break;
            }
            winnerItr=getLeafIter(winner);
            winnerItr.setCurrentMergedStart(mergedStart);
           // System.out.println("######"+Integer.toHexString(getLeafCode(winner))+":::"+winner+":::"+winnerItr.getCurrentInitStart()+":::"+winnerItr.getCurrentLen()+":::"+mergedStart);
            if(winnerCode!=getLeafCode(winner)) {
            	winnerCode=getLeafCode(winner);
            	finalIndexList.addCodeAndStartPos(winnerCode, mergedStart);
        		//System.out.println(Integer.toHexString(winnerCode)+":::"+mergedStart);
        	}
            mergedStart+=winnerItr.getCurrentLen();
        }
    	return finalIndexList;
    } 
}