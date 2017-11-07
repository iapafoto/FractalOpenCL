/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.tool;

import java.util.concurrent.Executor;
import javax.swing.SwingUtilities;

/**
 *
 * @author fguerry
 */
public class SwingExecutor implements Executor {

    private static final SwingExecutor instance = new SwingExecutor();
    
    private SwingExecutor() {
    }

    public static SwingExecutor getInstance() {
        return instance;
    }

    @Override
    public void execute(Runnable command) {
        if (SwingUtilities.isEventDispatchThread()) {
            command.run();
        } else {
            SwingUtilities.invokeLater(command);
        }
    }
    
}
