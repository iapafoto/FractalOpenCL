/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.tool;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.ListIterator;
import java.util.RandomAccess;
import java.util.Set;
import java.util.TreeSet;
import org.joda.time.DateTime;

/**
 *
 * @author sebastien
 * Recherches rapide dans des listes ordonnes
 * Recherche si l'element existe
 * Recherche l'element le plus proche
 * Recherche l'element precedent ou suivant
 * Recherche de l'element le plus proche dans un interval
 * Intersection rapide entre listes
 * etc.
 */
public abstract class SortedListTools {

    private static final int BINARYSEARCH_THRESHOLD = 5000;


    public static Integer findClosest(final List<Integer> list, final Integer key) {
        final Integer id = idOfClosest(list, key, 0, list.size()-1);
        if (id != null) {
            return list.get(id);
        }
        return null;
    }

    public static Integer idOfClosest(final List<Integer> list, final Integer key) {
        return idOfClosest(list, key, 0, list.size()-1);
    }
    public static Integer idOfClosest(final List<Integer> list, final Integer key, final int lowId, final int highId) {
        // Cas particuliers
        if (list == null || list.isEmpty() || lowId > highId) {
            return null;
        }
        // Evacuation des cas où la valeur serait hors de l'interval -----------
        if (list.get(lowId) >= key) {
            return lowId;
        }
        if (list.get(highId) <= key) {
            return highId;
        }

        int idx = binarySearch(list, key, lowId, highId);
        if (idx < 0) {
            idx = -(idx) - 1;
            if (idx != 0 && idx < list.size()) {
                return (key - list.get(idx - 1) <= list.get(idx) - key) ? idx - 1 : idx;
            }
            return null;
        }
        return idx;
    }

    public static Long findClosest(final List<Long> list, final Long key) {
        final Integer id = idOfClosest(list, key, 0, list.size()-1);
        if (id != null) {
            return list.get(id);
        }
        return null;
    }

    public static Long findClosest(final List<Long> list, final Long key, final Long start, final Long end) {
        if (list == null || list.isEmpty()) {
            return null;
        }
        // On recherche les id de depart et de fin
        // On restrein un peu l'interval
        final Integer idStart = idOfNextOrEqual(list, start);
        if (idStart == null) {
            return null;
        }
        final Integer idEnd = idOfPreviousOrEqual(list, end);
        if (idEnd == null) {
            return null;
        }

        final Integer id = idOfClosest(list, key, idStart, idEnd);
        if (id != null) {
            return list.get(id);
        }
        return null;
    }
    
    public static Integer idOfClosest(final List<Long> list, final Long key) {
        return idOfClosest(list, key, 0, list.size()-1);
    }
    public static Integer idOfClosest(final List<Long> list, final Long key, int lowId, int highId) {
        // Cas particuliers
        if (list == null || list.isEmpty() || lowId > highId) {
            return null;
        }
        // Evacuation des cas où la valeur serait hors de l'intervalle -----------
        if (list.get(lowId) >= key)  {
            return lowId;
        }
        if (list.get(highId) <= key) {
            return highId;
        }

      // Recherche par dychotomie ----------------------------------------------
        int idx = binarySearch(list, key, lowId, highId);
        if (idx < 0) {
            idx = -(idx) - 1;
            if (idx != 0 && idx < list.size()) {
                return (key - list.get(idx - 1) <= list.get(idx) - key) ? idx - 1 : idx;
            }
            return null;
        }
        return idx;
    }

    public static Double findClosest(final List<Double> list, final Double key) {
        if (list != null && key != null) {
            final Integer id = idOfClosest(list, key, 0, list.size()-1);
            if (id != null) {
                return list.get(id);
            }
        }
        return null;
    }
    public static Integer idOfClosest(final List<Double> list, final Double key) {
        return idOfClosest(list, key, 0, list.size()-1);
    }
    public static Integer idOfClosest(final List<Double> list, final Double key, final int low, final int high) {
        // Cas particuliers
        if (list == null || list.isEmpty() || low > high) {
            return null;
        }
      // Evacuation des cas où la valeur serait hors de l'intervalle -----------
        if (list.get(low) >= key) {
            return low;
        }
        if (list.get(high) <= key) {
            return high;
        }

        int idx = binarySearch(list, key, low, high);
        if (idx < 0) {
            idx = -(idx) - 1;
            if (idx != 0 && idx < list.size()) {
                return (key - list.get(idx - 1) <= list.get(idx) - key) ? idx - 1 : idx;
            }
            return null;
        }
        return idx;
    }

    public static DateTime findClosest(final List<DateTime> list, final DateTime key) {
        if (list != null && !list.isEmpty()) {
            Integer id = idOfClosest(list, key, 0, list.size()-1);
            if (id != null) {
                return list.get(id);
            }
        }
        return null;
    }
    public static Integer idOfClosest(final List<DateTime> list, final DateTime key) {
        if(list!=null && !list.isEmpty()){
            return idOfClosest(list, key, 0, list.size()-1);
        } else {
            return null;
        }
    }
    public static Integer idOfClosest(final List<DateTime> list, final DateTime key, final int low, final int high) {
        // Cas particuliers
        if (list == null || list.isEmpty() || key == null || low > high) {
            return null;
        }
        // Evacuation des cas où la valeur serait hors de l'intervalle -----------
        if (list.size() == 1) {
            return 0;
        }
        if (list.get(low).getMillis() >= key.getMillis()) {
            return low;
        }
        if (list.get(high).getMillis() <= key.getMillis()) {
            return high;
        }

        int idx = binarySearch(list, key, low, high);
        if (idx < 0) {
            idx = -(idx) - 1;
            if (idx != 0 && idx < list.size()) {
                return (key.getMillis()  - list.get(idx - 1).getMillis()  <= list.get(idx).getMillis() - key.getMillis() ) ? idx - 1 : idx;
            }
            return null;
        }
        return idx;
    }
    
    public static <T extends Comparable> T findNext(final List<? extends Comparable<? super T>> list, T key) {
        final Integer id = idOfNext(list,key);
        if (id != null) {
            return (T)list.get(id);
        }
        return null;
    }
    
    public static <T extends Comparable> T findNextOrEqual(final List<? extends Comparable<? super T>> list, final T key, T start, T end) {
        Integer id = idOfNextOrEqual(list,key, start, end);
        if (id != null) {
            return (T)list.get(id);
        }
        return null;
    }
    public static <T extends Comparable> Integer idOfNextOrEqual(final List<? extends Comparable<? super T>> list, T key, T start, T end) {
        // Evacuation des cas particuliers

        if (list == null || list.isEmpty() || // Pas de valeurs => Pas de solution
            key==null /*|| end==null || key.compareTo(end) > 0*/) { 
            return null;
        }
        if (end == null) {
            end = (T)(list.get(list.size()-1));
        }
        if (key.compareTo(end) > 0) {   // Apres la fin => Pas de solution
            return null;
        }
        if (key.compareTo(start) < 0) { // avant le debut, on ignore les solutions avant start => on recherche la valeur qui suit
            key = start;
        }
        Integer idx = idOfNextOrEqual(list, key);
        if (idx != null) {
            key = (T)list.get(idx);
            if (key.compareTo(end) <= 0) { // La condition de debut est garanti par les tests ci desssus
                return idx;
            }
        }
        return null;
    }
    
    /**
     * Retourne les nb valeurs precedent la valeur passe en parametre
     * @param <T>
     * @param list
     * @param key
     * @param nb
     * @return 
     */
    public static <T extends Comparable> List<T>  findPrevious(final List<? extends Comparable<? super T>> list, T key, int nb) {
        final Integer id = idOfPrevious(list,key);
        if (id != null) {
            final List<T> lst = new ArrayList<>();
            for (int i=0; i<nb; i++) {
                if (id-i<0) {
                    break;
                }    
                lst.add((T)list.get(id-i));
            }
            return lst;
        }
        return null;
    }
    
    /**
     * Retourne les nb valeurs suivant la valeur passe en parametre
     * @param <T>
     * @param list
     * @param key
     * @param nb
     * @return 
     */
    public static <T extends Comparable> List<T> findNext(final List<? extends Comparable<? super T>> list, T key, int nb) {
        final Integer id = idOfNext(list,key);
        if (id != null) {
            final List<T> lst = new ArrayList<>();
            for (int i=0; i<nb; i++) {
                if (id+i>=list.size()) {
                    break;
                }         
                lst.add((T)list.get(id+i));
            }
            return lst;
        }
        return null;
    }
    
    public static <T extends Comparable> T findPrevious(final List<? extends Comparable<? super T>> list, T key) {
        final Integer id = idOfPrevious(list,key);
        if (id != null) {
            return (T)list.get(id);
        }
        return null;
    }
    public static <T extends Comparable> T findPreviousOrEqual(final List<? extends Comparable<? super T>> list, T key, T start, T end) {
        final Integer id = idOfPreviousOrEqual(list,key, start, end);
        if (id != null) {
            return (T)list.get(id);
        }
        return null;
    }
    public static <T extends Comparable> Integer idOfPreviousOrEqual(final List<? extends Comparable<? super T>> list, T key, T start, T end) {
        // Evacuation des cas particuliers
        if (list == null || list.isEmpty() || // Pas de valeurs => Pas de solution
            key.compareTo(start) < 0) {  // Avant le debut => Pas de solution
            return null;
        }
        if (key.compareTo(end) > 0) { // Apres la fin, on ignore les solutions apres end => on recherche la valeur qui precede
            key = end;
        }
        Integer idx = idOfPreviousOrEqual(list, key);
        if (idx != null) {
            key = (T)list.get(idx);
            if (key.compareTo(start) >= 0) { // La condition de fin est garanti par les tests ci desssus
                return idx;
            }
        }
        return null;
    }

    public static <T> Integer idOfNextOrEqual(final List<? extends Comparable<? super T>> list, final T key) {
        if (list == null || list.isEmpty()) {
            return null;
        }
        if (list.get(0).compareTo(key) >= 0) {
            return 0;
        }
        int idx = Collections.binarySearch(list, key);
        if (idx < 0) {
            idx = -(idx) - 1;
            if (idx != 0 && idx < list.size()) {
                return idx;
            }
            return null;
        }
        return idx;
    }
    
    public static <T> Integer idOfNext(final List<? extends Comparable<? super T>> list, final T key) {
        if (list == null || list.isEmpty()) {
            return null;
        }
        if (list.get(0).compareTo(key) > 0) {
            return 0;
        }
        int idx = Collections.binarySearch(list, key);
        // Pas trouve on retourne la valeur suivante
        if (idx < 0) {
            idx = -(idx) - 1;
            if (idx != 0 && idx < list.size()) {
                return idx;
            }
            return null;
        }
        // trouve dans la liste, on retourne la valeur qui suit si elle existe
        if (idx < list.size()-1) {
            return idx+1;
        }
        return null;
    }
    
    public static <T> Integer idOfPreviousOrEqual(final List<? extends Comparable<? super T>> list, final T key) {
        if (list == null || list.isEmpty()) {
            return null;
        }
        if (list.get(list.size() - 1).compareTo(key) <= 0) {
            return list.size() - 1;
        }
        int idx = Collections.binarySearch(list, key);
        if (idx < 0) {
            idx = -(idx) - 1;
            if (idx != 0 && idx < list.size()) {
                return idx - 1;
            }
            return null;
        }
        return idx;
    }
    
    public static <T> Integer idOfPrevious(final List<? extends Comparable<? super T>> list, final T key) {
        if (list == null || list.isEmpty()) {
            return null;
        }
        // on cherche la position de la valeur
        int idx = Collections.binarySearch(list, key);
        // si negatif, la valeur n est pas dans la liste 
        if (idx < 0) {
            idx = -(idx) - 1; // id de la valeur suivante dans la liste
            if (idx != 0 && idx < list.size()) {
                return idx - 1;
            }
            return null;
        }
        // si > 0 alors on a trouve la valeur (=>il nous faut celle d'avant)
        if (idx > 0) {
            return (idx-1);
        } else {
            return null;
        }
    }

    public static <T> int binarySearch(final List<? extends Comparable<? super T>> list, final T key, int low, int high) {
        if (list instanceof RandomAccess || list.size()<BINARYSEARCH_THRESHOLD) {
            return indexedBinarySearch(list, key, low, high);
        }
        else {
            return iteratorBinarySearch(list, key, low, high);
        }
    }

    private static <T> int indexedBinarySearch(final List<? extends Comparable<? super T>> list, final T key, int low, int high) {
        int cmp, mid;
        while (low <= high) {
            mid = (low + high) >>> 1;
            cmp = list.get(mid).compareTo(key);

            if (cmp < 0) {
                low = mid + 1;
            } else if (cmp > 0) {
                high = mid - 1;
            } else {
                return mid;
            } // key found
        }
        return -(low + 1);  // key not found
    }

    private static <T> int iteratorBinarySearch(final List<? extends Comparable<? super T>> list, final T key, int low, int high)  {
        final ListIterator<? extends Comparable<? super T>> i = list.listIterator();
        int mid, cmp;
        while (low <= high) {
            mid = (low + high) >>> 1;
            cmp = get(i, mid).compareTo(key);

            if (cmp < 0) {
                low = mid + 1;
            } else if (cmp > 0) {
                high = mid - 1;
            } else {
                return mid;
            } // key found
        }
        return -(low + 1);  // key not found
    }

    private static <T> T get(final ListIterator<? extends T> i, final int index) {
        T obj = null;
        int pos = i.nextIndex();
        if (pos <= index) {
            do {
                obj = i.next();
            } while (pos++ < index);
        } else {
            do {
                obj = i.previous();
            } while (--pos > index);
        }
        return obj;
    }
    
    /**
     * Trouve l'intersection entre toutes les listes triees pass�es en parametre
     * @param listOList: liste de liste 
     * @return La listes des elements presentes dans toutes les listes
     */
    public static <T extends Comparable> List<T> intersection(List<T> list1, List<T> list2) {
        List<List<T>> lstLst = new ArrayList<>();
        lstLst.add(list1);
        lstLst.add(list2);
        return intersection(lstLst);
    }
//    public static <T extends Comparable> List<T> intersection(final Collection<List<T>> listOfList) {
//        return intersection(Lists.newArrayList(listOfList));
//    }
    
    public static <T extends Comparable> List<T> intersection(final List<List<T>> listOfList) {
        T candidate;

        // Definition d'un tableau de position dans les listes et initialisation au debut
        final int[] tabIndex = new int[listOfList.size()];
        final List<T> listCommon = new ArrayList<>();
        if (listOfList.isEmpty()) {
            return listCommon;
        } else if (listOfList.size() == 1) {
            return listOfList.get(0);
        }
        
        List<T> list;
        int lstId;
        
       // ajoute par isa pour macro:si une des listes est vide, on retourne une liste vide pour l'intersection
       for (List<T> lst : listOfList) {
            if (lst.isEmpty()) {
                return new ArrayList<>();
            }
        }
       
        while (true) { // On va forcement sortir car on ne fait qu'avancer dans les listes
            // Recuperation de la valeur de la premiere liste
            list = listOfList.get(0);
            if (tabIndex[0] >= list.size()) { // la liste est au bout => pas d'autre date possible
                return listCommon;
            }
            candidate = list.get(tabIndex[0]);
            lstId = 1;

            while (lstId < listOfList.size()) {
                list = listOfList.get(lstId);

                // On rattrape le retard sur cette liste (on pourrait le faire par dicotomie entre "posInList[lstId] et timeList.length" pour ameliorer encore l'efficacitee)
                while (list.get(tabIndex[lstId]).compareTo(candidate) < 0) {
                    tabIndex[lstId]++;
                    if (tabIndex[lstId] >= list.size()) { // Cette liste est au bout => pas d'autre valeur possibles
                        return listCommon;
                    }
                }

                // Si on a depass� la date => elle n'existe pas dans cette liste
                // => On definit une nouvelle date candidate avec la valeur suivante trouv�e
                // Et on relance l'algo � la premiere liste
                if (list.get(tabIndex[lstId]).compareTo(candidate) > 0) {
                    candidate = list.get(tabIndex[lstId]);
                    lstId = 0;
                    continue;
                }
                // Si on est arriv� jusque ici c'est que la date existe dans les lstId listes precedentes
                // => On va la rechercher dans la liste suivante
                lstId++;
            }
            // Si on arrive ici sans etre sorti avant c'est que l'on a une date valide pour toutes les listes
            listCommon.add(candidate);

            // On continu en repartant sur la prochaine date candidate de la premiere liste
            tabIndex[0]++;
        }

    }
    
    /**
     * Renvoi l'ensemble des valeurs qui sont dans l'une des listes 
     * Pas testé en performance
     * @param <T>
     * @param listOfList
     * @return 
     */
    public static <T extends Comparable> List<T> union(final Collection<List<T>> listOfList) {
        final Set<T> set = new TreeSet<>();
        for (List<T> lst : listOfList) {
            set.addAll(lst);
        }
        return new ArrayList<>(set); // TreeSet is ordered
    }
    
    public static <T extends Comparable> List<T> union(final List<T> lst1, final List<T> lst2) {
        final Set<T> set = new TreeSet<>(lst1);
        set.addAll(lst2);
        return new ArrayList<>(set); // TreeSet is ordered
    }
    /**
     * Retourne l'ensemble des dates de la deuxieme liste, les plus proches de celles de la 1ere liste 
     * TODO: possibilitee de le faire en (n+m) plutot qu'en m*log(n) en parcourant les deux liste en meme temps de maniere lineaire
     * @param master
     * @param lst
     * @return 
     */
    public static List<DateTime> closest(final List<DateTime> master, final List<DateTime> lst) {
        return closest(master, lst, null);
    }
    
    /**
     * Retourne l'ensemble des dates de la deuxieme liste, les plus proches de celles de la 1ere liste 
     * TODO: possibilitee de le faire en (n+m) plutot qu'en m*log(n) en parcourant les deux liste en meme temps de maniere lineaire
     * @param master
     * @param lst
     * @param maxInterval interval au dela duquel l'image est consideree comme invalide car trop eloignee
     * @return 
     */
    public static List<DateTime> closest(final List<DateTime> master, final List<DateTime> lst, final Long maxInterval) {
        if (master == null || lst == null) {
            return null;
        }
        final List<DateTime> closest = new ArrayList<>();
        if (master.isEmpty() || lst.isEmpty()) {
            return closest;
        }       
        if (master == lst) { // il y a de chances que ce soit la mme liste : on gagne au moins ca !
            return master;
        }
        int endId = lst.size()-1;
        int startId = 0;
        Integer id;
        // Restriction de la recherche a la partie contenant des valeurs interessantes
        id = idOfClosest(lst, master.get(master.size()-1), startId, endId);
        if (id != null) {
            endId = id;
        }
        for (DateTime value : master) {
            id = idOfClosest(lst, value, startId, endId);
            if (id != null) {
                if (maxInterval == null || Math.abs(value.getMillis() - lst.get(id).getMillis()) <= maxInterval) {
                    closest.add(lst.get(id));
                    startId = id+1; // Pour etre sur de ne plus reutiliser cet element là
                }
                if (startId > endId) {
                    break;
                }
            }
        }
        return closest;
    }

    public static boolean contains(final List<DateTime> synchronizedTimes, final DateTime time) {
        return (time != null) && time.equals(findClosest(synchronizedTimes, time));
    }

}
