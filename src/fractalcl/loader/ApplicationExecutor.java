/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.loader;


import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import net.sf.ehcache.util.NamedThreadFactory;

/**
 * 
 */
public class ApplicationExecutor {

    // === Propriétés ===
    /**
     * Instance du singleton
     */
    private static ApplicationExecutor instance = null;
    /**
     * Instance du RequestProcessor
     */
    private ListeningExecutorService service = null;
    private ScheduledExecutorService scheduledService = null;
    /**
     * durée de survie d'un thread excédentaire non utilisé en secondes (paramètre TimeUnit.SECONDS)..
     */
    private static final int keepAliveTime=1000;
    /**
     * factories de threads
     */
    private final ThreadFactory threadFactory;

    // === Méthodes publiques ===
    /**
     * Permet de récupérer l'instance du singleton
     * 
     * @return instance du singleton
     */
    public synchronized static ApplicationExecutor getInstance() {
        if (instance == null) {
            instance = new ApplicationExecutor();
        }
        return instance;
    }
    
    // === Méthodes privées ou protégées ===
    /**
     * Constructeur privé
     */
    private final BlockingQueue<Runnable> workQueue;
    private final ThreadPoolExecutor threadExecutor;
    private final int maxSize;
    
    private ApplicationExecutor() { 
        threadFactory = new NamedThreadFactory("synopsisThreads");
        workQueue = new LinkedBlockingQueue<>(3000);
        maxSize = 40;
        threadExecutor = new ThreadPoolExecutor(20, maxSize, keepAliveTime, TimeUnit.SECONDS, workQueue,threadFactory);
        service = MoreExecutors.listeningDecorator(threadExecutor);
        scheduledService = Executors.newScheduledThreadPool(20);
    }

    // === Méthodes triviales ou autogénérées ===
    public ListeningExecutorService getExecutor() {
        return service;
    }
    public ScheduledExecutorService getScheduledExecutor(){
        return scheduledService;
    }
    public List<String> getInfos(){
        List<String> infos= new ArrayList<>();
        infos.add("-----Synopsis Thread pool----");
        infos.add("config: max pool size : "+threadExecutor.getMaximumPoolSize());
        infos.add("config: core pool size : "+threadExecutor.getCorePoolSize());
        infos.add("config: queue capacity (from configuration) : "+maxSize);
        infos.add("pool size : "+threadExecutor.getPoolSize());
        infos.add("nb active tasks : "+threadExecutor.getActiveCount());
        infos.add("max nb tasks : "+threadExecutor.getLargestPoolSize());
        infos.add("max ever nb tasks : "+threadExecutor.getLargestPoolSize());
        infos.add("nb total executed tasks : "+threadExecutor.getTaskCount());
        infos.add("queue size : "+threadExecutor.getQueue().size());
        infos.add("queue remaining capacity : "+threadExecutor.getQueue().remainingCapacity());
        return infos;
    }
    /*
    public ApplicationExecutorsInfos getExecutorStats(){
        return new  ApplicationExecutorsInfos(threadExecutor.getMaximumPoolSize(),
                 threadExecutor.getCorePoolSize(),
                 maxSize,
                threadExecutor.getPoolSize(),
                threadExecutor.getActiveCount(),
                threadExecutor.getLargestPoolSize(),
                threadExecutor.getTaskCount(),
                threadExecutor.getQueue().size());
    }*/
}
