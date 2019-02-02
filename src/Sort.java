import java.util.Arrays;

/*
Sorting Algorithms
- Bubble Sort
- Selection Sort
- Insertion Sort
- Heap Sort
- Merge Sort
- Quick Sort
- Quick Sort (3-Way Partition)
*/

public class Sort {

    public static void bubbleSort(int array[]) {
        if (array == null || array.length == 0)
            return;

        // The largest element is in correct sorted order after each pass
        for (int i = (array.length - 1); i >= 0; i--) {
            for (int j = 1; j <= i; j++) {
                if (array[j - 1] > array[j]) {
                    int temp = array[j - 1];
                    array[j - 1] = array[j];
                    array[j] = temp;
                }
            }
        }
    }

    public static void selectionSort(int[] array) {
        if (array == null || array.length == 0)
            return;

        for (int i = 0; i < array.length - 1; i++) {
            // Linear search for the minimum
            int min = i;
            for (int j = i + 1; j < array.length; j++) {
                if (array[j] < array[min]) {
                    min = j;
                }
            }

            // Put minimum into correct place
            if (min != i) {
                int temp = array[i];
                array[i] = array[min];
                array[min] = temp;
            }
        }
    }

    public static void insertionSort(int[] array) {
        if (array == null || array.length == 0)
            return;

        for (int i = 1; i < array.length; i++) {
            // Take an unsorted element and insert it into correct sorted spot
            for (int j = i; j > 0 && array[j] < array[j - 1]; j--) {
                int temp = array[j];
                array[j] = array[j - 1];
                array[j - 1] = temp;
            }
        }
    }

    private static void heapSort(int[] array) {
        Heapify(array);

        // The array from [0, ... , i] is the heap
        // The array from [i+1, ... , end] is in sorted order
        for (int i = array.length - 1; i > 0; i--) {
            int temp = array[0];
            array[0] = array[i];
            array[i] = temp;

            siftDown(array, 0, i - 1);
        }
    }

    private static void Heapify(int[] array) {
        // Root Node is at index 0
        // Nth Node is at index (n + 1)
        // Left Child is at index (2n + 1)
        // Right Child is at index (2n + 2)

        // Start with the parent of the last element
        // Sift down the node to its proper place
        for (int node = (array.length - 2) / 2; node >= 0; node--) {
            siftDown(array, node, array.length - 1);
        }
    }

    private static void siftDown(int[] array, int root, int end) {
        // Continually swap the root with the greatest child until there's no
        // more children or the root is greater than the children
        while (2 * root + 1 <= end) {
            int swap = root;
            int child = 2 * root + 1; // left child

            // Check left child
            if (array[child] > array[swap]) {
                swap = child;
            }

            // Check right child
            if ((child + 1 <= end) && (array[child + 1] > array[swap])) {
                swap = child + 1;
            }

            if (swap == root) {
                // The root is greater than its children
                return;
            } else {
                int temp = array[root];
                array[root] = array[swap];
                array[swap] = temp;

                root = swap;
            }
        }
    }

    public static void mergeSort(int[] array) {
        if (array.length <= 1)
            return;

        int[] left = new int[array.length / 2];
        int[] right = new int[array.length - left.length];

        for (int i = 0; i < left.length; i++) {
            left[i] = array[i];
        }

        for (int i = 0; i < right.length; i++) {
            right[i] = array[left.length + i];
        }

        mergeSort(left);
        mergeSort(right);
        merge(left, right, array);
    }

    private static void merge(int[] first, int[] second, int[] array) {
        int indexFirst = 0;
        int indexSecond = 0;
        int j = 0;

        while (indexFirst < first.length && indexSecond < second.length) {
            if (first[indexFirst] < second[indexSecond]) {
                array[j] = first[indexFirst];
                indexFirst++;
            } else {
                array[j] = second[indexSecond];
                indexSecond++;
            }

            j++;
        }

        while (indexFirst < first.length) {
            array[j] = first[indexFirst];
            indexFirst++;
            j++;
        }

        while (indexSecond < second.length) {
            array[j] = second[indexSecond];
            indexSecond++;
            j++;
        }
    }

    public static void quickSort(int[] array) {
        quickSort(array, 0, array.length - 1);
    }

    private static void quickSort(int[] array, int lower, int higher) {
        if (lower < higher) {
            int pivotIndex = partition(array, lower, higher);

            quickSort(array, lower, pivotIndex - 1);
            quickSort(array, pivotIndex + 1, higher);
        }
    }

    private static int partition(int[] array, int lower, int higher) {
        // Pick a pivot between lower index and higher index
        int pivotIndex = (int)(Math.random() * (higher - lower))
                + lower;

        // Place pivot at the end of the subarray
        int temp = array[pivotIndex];
        array[pivotIndex] = array[higher];
        array[higher] = temp;

        pivotIndex = higher;
        int pivot = array[pivotIndex];

        // Partitioning
        int i = lower;
        for (int j = lower; j < higher; j++) {
            if (array[j] <= pivot) {
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;

                i++;
            }
        }

        // Move the pivot into position i
        temp = array[i];
        array[i] = array[pivotIndex];
        array[pivotIndex] = temp;

        return i;
    }

    public static void quickSort3(int[] array) {
        // Quick Sort with a 3-Way Partition: [ < Pivot | = Pivot | > Pivot ]
        quickSort3(array, 0, array.length - 1);
    }

    public static void quickSort3(int[] array, int lower, int higher) {
        if (lower < higher) {
            // Pick a pivot between lower index and higher index
            int pivotIndex = (int)(Math.random() * (higher - lower))
                    + lower;

            // Place pivot at the start of the (sub)array
            int temp = array[pivotIndex];
            array[pivotIndex] = array[lower];
            array[lower] = temp;

            pivotIndex = lower;
            int pivot = array[pivotIndex];

            // The upper index of the "less than pivot" partition
            int lt = lower;
            // The lower index of the "greater than pivot" partition
            int gt = higher;

            // Partitioning
            int i = lower + 1;
            while (i <= gt) {
                if (array[i] < pivot) {
                    temp = array[lt];
                    array[lt] = array[i];
                    array[i] = temp;

                    i++;
                    lt++;

                } else if (array[i] > pivot) {
                    temp = array[gt];
                    array[gt] = array[i];
                    array[i] = temp;

                    gt--;

                } else {
                    i++;
                }
            }

            // Recursively call on the left and right partitions
            quickSort3(array, lower, lt - 1);
            quickSort3(array, gt + 1, higher);
        }
    }

    public static boolean isSorted(int[] array) {
        for (int i = 0; i < array.length - 1; i++) {
            if (array[i] > array[i + 1]) {
                return false;
            }
        }
        return true;
    }

    // Initialize array to integers between 0 and n
    public static void randomArray(int[] array, int n) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (int) (Math.random() * n);
        }
    }

    public static void main(String[] args) {
        int[] array = new int[10000];
        int n = 9000;
        long startTime, endTime;

        // Bubble Sort
        randomArray(array, n);
        startTime = System.nanoTime();
        bubbleSort(array);
        endTime = System.nanoTime();
        System.out.println("Bubble Sort:");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);

        // Selection Sort
        randomArray(array, n);
        startTime = System.nanoTime();
        selectionSort(array);
        endTime = System.nanoTime();
        System.out.println("\nSelection Sort:");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);

        // Insertion Sort
        randomArray(array, n);
        startTime = System.nanoTime();
        insertionSort(array);
        endTime = System.nanoTime();
        System.out.println("\nInsertion Sort:");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);

        // Heap Sort
        randomArray(array, n);
        startTime = System.nanoTime();
        heapSort(array);
        endTime = System.nanoTime();
        System.out.println("\nHeap Sort:");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);

        // Merge Sort
        randomArray(array, n);
        startTime = System.nanoTime();
        mergeSort(array);
        endTime = System.nanoTime();
        System.out.println("\nMerge Sort:");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);

        // Quick Sort
        randomArray(array, n);
        startTime = System.nanoTime();
        quickSort(array);
        endTime = System.nanoTime();
        System.out.println("\nQuick Sort:");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);

        // Quick Sort (3 way partition)
        randomArray(array, n);
        startTime = System.nanoTime();
        quickSort3(array);
        endTime = System.nanoTime();
        System.out.println("\nQuick Sort (3-Way Partition):");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);

        // Quick Sort (3 way partition)
        randomArray(array, n);
        startTime = System.nanoTime();
        Arrays.sort(array);
        endTime = System.nanoTime();
        System.out.println("\nJava Arrays.sort:");
        System.out.println("(" + isSorted(array) + "): "
                + (endTime - startTime) / 1000000000.0);
    }
}