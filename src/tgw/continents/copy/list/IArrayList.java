package tgw.continents.copy.list;

import it.unimi.dsi.fastutil.ints.IntArrayList;

public class IArrayList extends IntArrayList {

    @Override
    public int set(int index, int k) {
        if (index > this.size) {
            throw new IndexOutOfBoundsException("Index (" + index + ") is greater than or equal to list size (" + this.size + ")");
        }
        int old = this.a[index];
        this.a[index] = k;
        return old;
    }
}
