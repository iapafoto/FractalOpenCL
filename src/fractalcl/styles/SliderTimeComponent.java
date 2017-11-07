package fractalcl.styles;



import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.List;
import java.util.Set;
import fractalcl.styles.FlyingToolTip.FlyingToolTipListener;
import fractalcl.tool.SortedListTools;
import org.joda.time.DateTime;
import org.joda.time.DateTimeZone;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;

/**
 *
 * @author sebastien.durand
 */
public final class SliderTimeComponent extends SliderComponentForTimeSlider implements FlyingToolTipListener, MouseMotionListener, MouseListener, FocusListener, MouseWheelListener, KeyListener {
    
    enum RefreshMode {
        NONE, SLIDING, STATIC
    }
    
    public abstract class TimeConstant {

        public static final long SECOND = 1000;
        public static final long MINUTE = 60 * SECOND;
        public static final long HOUR = 60 * MINUTE;
        public static final long DAY = 24 * HOUR;
        public static final long WEEK = 7 * DAY;
        // Valeures speciales (dans ce cas la valeur correspond au facteur d'echelle relative)
        public static final long MONTH = 31 * DAY;
        public static final long YEAR = 366 * DAY;
    }
    
    /**
     * Evenements
     */
    public static final String 
            EVT_LIMITS_CHANGED  = "LIMITS_CHANGED",
            EVT_LIMITS_CHANGING = "LIMITS_CHANGING";
    /**
     * Dessins des curseurs
     */
    public final static Shape 
            SHAPE_LIMIT_MIN  = new Polygon(new int[] {0, 0,-7}, new int[] {0,-14, 0}, 3),
            SHAPE_LIMIT_MAX  = new Polygon(new int[] {0, 0, 7}, new int[] {0,-14, 0}, 3),
            SHAPE_LIMIT_AUTO = new Polygon(new int[] {-2,-2, 0, 5, 0}, new int[] { 1,-13,-13,-6, 1}, 5);
    
    /**
     * Couleurs utilisés pour le dessin
     */
    public final static Color
            COLOR_SLIDER = new Color(64,64,64),
            COLOR_VALID_FRAME = new Color(0, 255, 255, 128),
            COLOR_FILTERED_FRAME = new Color(255, 100, 50),
            COLOR_CURSOR_LINE = new Color(255,255,255,32),
            COLOR_AVAILABLE_FRAME = ColorTools.alpha(Color.lightGray, 64),
            COLOR_ALPHA_WHITE_20 = ColorTools.alpha(Color.white, 20),
            COLOR_ALPHA_WHITE_40 = ColorTools.alpha(Color.white, 40),
            COLOR_ALPHA_WHITE_60 = ColorTools.alpha(Color.white, 60),
            COLOR_ALPHA_WHITE_80 = ColorTools.alpha(Color.white, 80),    
            COLOR_PAST = new Color(0x18ed4a),
            COLOR_FUTURE = new Color(0xffe200);
        
    /**
     * Formatter de date pour les tooltips et les popups
     */
    public static DateTimeFormatter dateFormatter = DateTimeFormat.forPattern("EEE dd/MM HH:mm").withZoneUTC();

    /**
     * Mode de refresh de l'animation associee
     */
    private RefreshMode refreshMode = RefreshMode.NONE;
    
    /**
     * La liste des dates valides dans le slider
     */
    private List<DateTime> 
            availableDates = null, // Liste des dates disponibles
            validDates     = null, // Liste des dates selectionnee dans l'animation
            validTimesWithoutFiltredFrames = null, // Liste des dates affichees dans l'animation
            inCacheDates   = null; // Liste des dates en cache
    private Set<DateTime> specialDates = null; // Liste des dates filtrées
    private Double 
            dateMinAsValue, // Borne de debut d'animation en UNIXTIME
            dateMaxAsValue; // Borne de fin d'animation en UNIXTIME
    
    /**
     * Flags internes pour gerer les selection possible en fonction des modes
     */
    private boolean 
            isMinLimitFixed  = false, 
            isMaxLimitFixed  = false,
            minLimitDragging = false,
            maxLimitDragging = false,
            minLimitSelected = false,
            maxLimitSelected = false,
            sliderSelected   = false,
            minLimitOver     = false,
            maxLimitOver     = false,
            isZoomedOnValid  = false;
   

    
    /**
     * Popup-Tooltip indiquant la date selectionnee à un moment donné
     */
    private FlyingToolTip timePopup;
    
    
    public SliderTimeComponent() {
        this(new DateTime(DateTimeZone.UTC), new DateTime(DateTimeZone.UTC).plus(TimeConstant.DAY), new DateTime(DateTimeZone.UTC).plus(TimeConstant.DAY/2));
    }

    public SliderTimeComponent(final DateTime limitMin, final DateTime limitMax, final DateTime value) {
        this(limitMin, limitMax);
        setValue(value.getMillis());
    }

    /**
     * Constructeur complet
     * @param limitMin
     * @param limitMax 
     */
    public SliderTimeComponent(final DateTime limitMin, final DateTime limitMax) {
        super(0, limitMin.getMillis(), limitMax.getMillis(), Color.white, Color.white);

        setMaximumSize(new Dimension(2000, 20));
        setPreferredSize(new Dimension(2000, 20));
        setMinimumSize(new Dimension(60, 20));
        
        timePopup = new FlyingToolTip(this);
        
        addMouseWheelListener(this);
        addFocusListener(this);
        addKeyListener(this);
    }

    /**
     * Changement des bornes de l'animation
     * @param limitMin
     * @param limitMax 
     */
    public void setMinAndMax(final DateTime limitMin, final DateTime limitMax) {
        if (!minLimitDragging && !minLimitSelected) {
            this.dateMinAsValue = (limitMin != null) ? (double)limitMin.getMillis() : null;
        }
        if (!maxLimitDragging && !maxLimitSelected) {
            this.dateMaxAsValue = (limitMax != null) ? (double)limitMax.getMillis() : null;
        }
    }

    /**
     * Changement du mode de refresh de l'animation
     * @param mode 
     */
  //  public void setRefreshMode(final RefreshMode mode) {
  //      refreshMode = mode;
  //      isMaxLimitFixed = (/*isZoomedOnValid ||*/ (refreshMode != RefreshMode.NONE));
  //      isMinLimitFixed = (/*isZoomedOnValid ||*/ (refreshMode == RefreshMode.STATIC));
  //  }
    
    /**
     * Changement du formatter de date 
     * Pourrait etre utilisé pour changer de fuseau horaire par exemple
     * @param dateFormatter 
     */
    public void setDateFormatter(final DateTimeFormatter dateFormatter) {
        this.dateFormatter = dateFormatter;
    }

    /**
     * Changement des dates utilisees dans la barre de temps
     * @param availableDates liste des dates existantes
     * @param validDates liste des dates utilisee par l'animation
     * @param special liste des dates filtrees
     * @param inCache liste des dates qui sont dans le cache
     */
    public void setValidDates(final List<DateTime> availableDates, final List<DateTime> validDates, final List<DateTime> validTimesWithoutFiltredFrames,
                              final Set<DateTime> special, final List<DateTime> inCache) {
        this.availableDates = availableDates;
        this.validDates = validDates;
        this.validTimesWithoutFiltredFrames = validTimesWithoutFiltredFrames;
        this.specialDates = special;
        this.inCacheDates = inCache;
        updateDisplayLimit();
        repaint();
    }

    public void addDatesToCache(final List<DateTime> newInCache) {
        // On ajoute dans l'ordre en evitant les doublons
        inCacheDates = SortedListTools.union(inCacheDates, newInCache);
        repaint();
    }
    
    /**
     * Changement de la date selectionné dans la barre de temps
     * @param value 
     */
    public void setDate(final DateTime value) {
        super.setValue(value != null ? value.getMillis() : 0);
        repaint();
    }

    /**
     * Recuperation de la date selectionnee (forcement dans la liste de dates valid)
     * @return 
     */
    public DateTime getDate() {
        return SortedListTools.findClosest(validDates, new DateTime((long)getValue()));
    }

    /**
     * Transformation d'une valeur du slider en texte
     * On choisit, parmis les dates valides, la date la plus proche de la valeur passée
     * Rq: Utilisé par la classe mère
     * @param value
     * @return 
     */
    @Override
    protected String valueAsString(final double value) {
        final DateTime date = SortedListTools.findClosest(validDates, new DateTime((long)value));
        return dateAsString(date);
    }
    
    /**
     * Transformation d'une valeur "libre" du slider en texte
     * On choisit, parmis toutes les dates existantes, la date la plus proche de la valeur passée
     * @param value
     * @return 
     */
    protected String valueOfAvailableAsString(final double value) {
        final DateTime date = SortedListTools.findClosest(availableDates, new DateTime((long)value));
        return dateAsString(date);
    }
    
    /**
     * Formattage des dates en texte à afficher dans l'IHM
     * @param date
     * @return 
     */
    public static String dateAsString(final DateTime date) {
        return (date != null) ? (dateFormatter != null) ? dateFormatter.print(date) : date.toString() : null;
    }
    
    /**
     * Transformation d'une position en une date (UNIXTIME) valide de la liste
     * @param pos
     * @return 
     */
    @Override
    protected double positionToValue(final int pos) {
        double result = super.positionToValue(pos);
        if (validDates != null) { // On recherhce une date parmis celles de la liste
            DateTime d = new DateTime((long) result);
            d = SortedListTools.findClosest(validDates, d);
            if (d != null) {
                return d.getMillis();
            }
        }
        return result;
    }

    /**
     * Recuperation de la valeur a envoye au propertyChangeListener lors des modification
     * Utilisé par la classe mère
     * @return 
     */
    @Override
    protected Object getValueForEvent() {
        return getDate();
    }
    
    /**
     * Dessin de tout le composant excepté le curseur de date
     * @param g
     */
    @Override
    protected void drawSlider(final Graphics2D g) {
        colorLimitMin = colorLimitMax = Color.lightGray;
        Graphics2D g2 = (Graphics2D) g.create();
        // Avoir de l'anti-aliasing sur l'intégralité de la méthode accélère les traitements (surprenant mais bon, 
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setColor(COLOR_SLIDER);
        fillSlider(g2);
        final int sliderY = sliderY();
        int posCurrentTime = Integer.MAX_VALUE;

        final int y1DrawLine = sliderY + 1;
        final int y2DrawLine = sliderY + SLIDER_H - 1;

        g2.setColor(COLOR_AVAILABLE_FRAME);
        if (availableDates != null) {
            // Affichage des donnees potentielement disponibles
            int lastX = Integer.MIN_VALUE;
            for (final DateTime d : availableDates) {
                final int posX = this.valueToPosition(d.getMillis());
                if (lastX != posX) {
                    g2.drawLine(posX, y1DrawLine, posX, y2DrawLine);
                }
                lastX = posX;
            }

            // Affichage de la date courante si les dates dispo incluent la date courante
            final long time = System.currentTimeMillis();
            final int s = availableDates.size();
            if (s > 1 && availableDates.get(0).getMillis() < time &&
                         availableDates.get(s-1).getMillis() > time) {
                posCurrentTime = valueToPosition(time);
                g2.setColor(Color.orange);
                g2.fillRect(posCurrentTime, sliderY-3, 2, SLIDER_H + 7);
            }
        }
   
        if (validDates != null) {
            if (validDates.isEmpty()) {
                g2.setColor(Color.white);
                fillSlider(g2);
            } else {
                g2.setColor(COLOR_VALID_FRAME);
                int lastX = Integer.MIN_VALUE;
                for (final DateTime d : validDates) {
                    final int posX = this.valueToPosition(d.getMillis());
                    if (posX != lastX) {
                        g2.drawLine(posX, y1DrawLine, posX, y2DrawLine);
                    }
                    lastX = posX;
                }
                g2.setColor(COLOR_PAST);
                lastX = Integer.MIN_VALUE;
                for (final DateTime d : inCacheDates) {
                    final int posX = this.valueToPosition(d.getMillis());
                    if (posX != lastX) {
                        if (posX > posCurrentTime) {
                            g2.setColor(COLOR_FUTURE);
                            posCurrentTime = Integer.MAX_VALUE; // pour ne pas le refaire a chaque fois
                        }
                        g2.drawLine(posX, y1DrawLine, posX, y2DrawLine);
                    }
                    lastX = posX;
                }

                if (specialDates != null) {
                    g2.setColor(COLOR_FILTERED_FRAME);//,Icons.getButtonColor());
                    lastX = Integer.MIN_VALUE;
                    for (final DateTime d : specialDates) {
                        final int posX = this.valueToPosition(d.getMillis());
                        if (posX != lastX) {
                            g2.drawLine(posX, y1DrawLine, posX, y2DrawLine);
                        }
                        lastX = posX;
                    }
                }
            }
        }
        
        if (isZoomedOnValid) {
            g2.setColor(COLOR_CURSOR_LINE);
            int w = getWidth();
            g2.drawLine(CURSOR_W, sliderY, w-CURSOR_W, sliderY);
            g2.drawLine(CURSOR_W, sliderY+SLIDER_H, w-CURSOR_W, sliderY+SLIDER_H);
        }
        


                
        g2.setColor(COLOR_ALPHA_WHITE_40);
        g2.fillRect(borderLeftRight, sliderY+1, sliderWidth(), 1);
        g2.setColor(COLOR_ALPHA_WHITE_60);
        g2.fillRect(borderLeftRight, sliderY+2, sliderWidth(), 1);
        g2.setColor(COLOR_ALPHA_WHITE_20);
        g2.fillRect(borderLeftRight, sliderY+3, sliderWidth(), 1);
        g2.setColor(COLOR_ALPHA_WHITE_80);
        g2.fillRect(borderLeftRight, sliderY+SLIDER_H-1, sliderWidth(), 1);
        g2.setColor(Color.black);
        g2.drawRoundRect(borderLeftRight, sliderY, sliderWidth(), SLIDER_H, ARC_SIZE, ARC_SIZE);
        
        // Dessin des limites de validitee min et max
        Color color;
        Shape shape, shapeCursor;

        if (refreshMode != RefreshMode.STATIC) {
            if (dateMinAsValue != null) {
                if (!isMinLimitFixed) {
                    color = (minLimitOver || minLimitDragging || minLimitSelected) ? Color.red : Color.orange;
                } else {
                    color = Color.gray;
                } 
                g2.setColor(color);
                final int posX = valueToPosition(this.dateMinAsValue);
                shapeCursor = (refreshMode == RefreshMode.SLIDING) ? SHAPE_LIMIT_AUTO : SHAPE_LIMIT_MIN;
                shape = AffineTransform.getTranslateInstance(posX+1, sliderY+9).createTransformedShape(shapeCursor);
                g2.fill(shape);
                g2.setColor(color.darker().darker().darker());
                // Parceque java decale le fill par rapport au draw !!!
               // shape = AffineTransform.getTranslateInstance(posX, sliderY+9).createTransformedShape(shapeCursor);
                g2.draw(shape);
            }

            if (dateMaxAsValue != null) {
                if (!isMaxLimitFixed) {
                    color = (maxLimitOver || maxLimitDragging || maxLimitSelected) ? Color.red : Color.orange;
                } else {
                    color = Color.gray;
                }
                final int posX = valueToPosition(dateMaxAsValue);
                shapeCursor = (refreshMode == RefreshMode.NONE) ? SHAPE_LIMIT_MAX : SHAPE_LIMIT_AUTO;
                shape = AffineTransform.getTranslateInstance(posX, sliderY+9).createTransformedShape(shapeCursor);
                g2.setColor(color);
                g2.fill(shape);
                g2.setColor(color.darker().darker().darker());
                g2.draw(shape);
            }
        } else { // RefreshMode.STATIC on affiche que la derniere borne en dedoublee
            if (dateMaxAsValue != null) {
                color = Color.yellow;
                final int posX = valueToPosition(dateMaxAsValue);
                shapeCursor = SHAPE_LIMIT_AUTO;
                
                shape = AffineTransform.getTranslateInstance(posX+2, sliderY+9).createTransformedShape(shapeCursor);
                g2.setColor(color);
                g2.fill(shape);
                g2.setColor(color.darker().darker().darker());
                g2.draw(shape);

                shape = AffineTransform.getTranslateInstance(posX-1, sliderY+9).createTransformedShape(shapeCursor);
                g2.setColor(color);
                g2.fill(shape);
                g2.setColor(color.darker().darker().darker());
                g2.draw(shape);
            }  
        }
        g2.dispose();
    }
    
    /**
     * Dessin du curseur de date
     * @param g2
     * @param cursorX 
     */
    @Override
    protected void drawCursor(final Graphics2D g2, int cursorX) {
        // En mode refresh statique : pas de curseur, on est figé sur la derniere date  
  //      if (!RefreshMode.STATIC.equals(refreshMode)) {
            final Color fillColor;
            if (specialDates != null && specialDates.contains(getDate())) {
                fillColor = Color.red;
            } else {
                fillColor = (cursorX > valueToPosition(System.currentTimeMillis())) ? COLOR_FUTURE : COLOR_PAST;             
            }
            drawCursor(g2, fillColor, Color.black, cursorX);
    //    }
    }
    
    /**
     * Donne la couleur associe au slider suivant sa date
     * @param date
     * @return 
     */
    public static Color dateAsColor(final long date) {
        return (date > System.currentTimeMillis()) ? COLOR_FUTURE : COLOR_PAST; 
    }

    
    /**
     * Regroupement du code pour mouvement du slider limit min
     * Met a jour la popup et le dessin du slider
     * @param date
     * @param isChanging indique si la veuleur est en cour de changement
     */
    private void moveMinLimit(final DateTime date, final boolean isChanging) {
         if (date != null) {
            timePopup.move(valueToPosition(date.getMillis()), dateAsString(date), null);
            onMinMaxChanged(date.getMillis(), dateMaxAsValue, isChanging);
            repaint();
        }
    }
    
    /**
     * Regroupement du code pour mouvement du slider limit max
     * Met a jour la popup et le dessin du slider
     * @param date
     * @param isChanging indique si la veuleur est en cour de changement
     */
    private void moveMaxLimit(final DateTime date, final boolean isChanging) {
         if (date != null) {
            timePopup.move(valueToPosition(date.getMillis()), dateAsString(date), null);
            onMinMaxChanged(dateMinAsValue, date.getMillis(), isChanging);
            repaint();
        }
    }
    
    /**
     * Regroupement du code pour mouvement du slider de temps
     * Met a jour la popup et le dessin du slider
     * @param date
     * @param isChanging indique si la valeur est en cour de changement
     */
    private void moveSlider(final DateTime date, final boolean isChanging) {
        if (date != null) {
            super.setValue(date.getMillis());
            if (isChanging) {
                timePopup.move(valueToPosition(value), valueAsString(value), dateAsColor((long)value));
                onDragging();
            } else {
                firePropertyChange(VALUE_CHANGED, null, date);
            }
            repaint();
            timePopup.move(valueToPosition(value), valueAsString(value), dateAsColor((long)value));
        }
    }
    
    /**
     * Appele apres mouvement des slider min ou max
     * Corrige les limites puis communique les nouvelles valeurs aux propertyChangeListeners
     * @param min
     * @param max
     * @param isChanging 
     */
    private void onMinMaxChanged(double min, double max, final boolean isChanging) {
     //   min = validateValueInLimits(min);
     //   max = validateValueInLimits(max);
        dateMinAsValue = Math.min(min, max);
        dateMaxAsValue = Math.max(min, max);
        firePropertyChange(isChanging ? EVT_LIMITS_CHANGING : EVT_LIMITS_CHANGED, null, new DateTime[] {
                new DateTime((long)(double)dateMinAsValue),
                new DateTime((long)(double)dateMaxAsValue)});
    }

    
    /**
     * Mise a jour des limites 
     * On se base sur les valeurs valides min et max pour plus declartee
     */
    private void updateDisplayLimit() {
        if (isZoomedOnValid) {
            if (validDates != null && !validDates.isEmpty()) {
                super.setLimits(
                        validDates.get(0).getMillis(),
                        validDates.get(validDates.size()-1).getMillis());
            }
        } else {
            if (availableDates != null && !availableDates.isEmpty()) {
                super.setLimits(
                        availableDates.get(0).getMillis(),
                        availableDates.get(availableDates.size()-1).getMillis());
            }
        }
    }
    
    /**
     * Recuperation du texte de tooltip
     * Quand la poput est affichee, on le masque afin de ne pas interferer avec elle
     * @param e
     * @return 
     */
    @Override
    public String getToolTipText(MouseEvent e) {
        if (minLimitSelected || maxLimitSelected || sliderSelected) {
            return null;
        }
        // #SYN-2945 (Jira) // Pas de jauge en dehors de l'interval
        final double v = super.positionToValue(e.getX());
        if (dateMinAsValue==null || dateMaxAsValue==null || v<dateMinAsValue || v>dateMaxAsValue) {
            return null;
        }
        return valueAsString(this.positionToValue(e.getX()));  // selection parmis les valeurs selectionables
    //  return valueOfAvailableAsString(super.positionToValue(e.getX())); // selection parmis toutes les valeurs possibles
    }
    
    /**
     * Appelé quand le popup se ferme à sa propre initiative
     * Exemple : sur mouvement ou resize de la fenetre mére
     */
    @Override
    public void onPopupHide() {
        minLimitSelected = maxLimitSelected = sliderSelected = false;
        repaint();
    }
    
    /**
     * Active / Desactive le zoom sur la partie valide des données
     * @param zoomed 
     */
    public void setZoomOnValid(final boolean zoomed) {
        isZoomedOnValid = zoomed;
        if (isZoomedOnValid) {
            borderLeftRight = 20;
        //    isMinLimitFixed = true;
        //    isMaxLimitFixed = true;
        } else {
            borderLeftRight = CURSOR_W; 
        //    isMinLimitFixed = false;
            isMaxLimitFixed = !(refreshMode == RefreshMode.NONE);
        }
        updateDisplayLimit();
        repaint();
    }
    
    
//==============================================================================    
//          +----------------------------------------+
//          |       GESTION DE LA SOURIS             |        
//          +----------------------------------------+

    /**
     * Sur mouse pressed on selectionne l'un des sliders pour ensuite le dragger
     * @param e 
     */
    @Override
    public void mousePressed(MouseEvent e) {
        this.requestFocus();
        
        minLimitSelected = maxLimitSelected = sliderSelected = false;
        
        // Aucune selection possible en mode veille statique
        if (RefreshMode.STATIC.equals(refreshMode)) {
            return;
        }
        
        final Integer 
                posMinLimit = (dateMinAsValue != null) ? valueToPosition(dateMinAsValue) : null,
                posMaxLimit = (dateMaxAsValue != null) ? valueToPosition(dateMaxAsValue) : null;

        if (!isMinLimitFixed && posMinLimit != null && e.getX() < posMinLimit+5 && e.getX() > posMinLimit - 10) {
            minLimitDragging = true;
            timePopup.show(e.getX(), valueOfAvailableAsString(super.positionToValue(e.getX())), Color.red);
            
        } else if (!isMaxLimitFixed && posMaxLimit != null && e.getX() > posMaxLimit-5 && e.getX() < posMaxLimit + 10) {
            maxLimitDragging = true;
            timePopup.show(e.getX(), valueOfAvailableAsString(super.positionToValue(e.getX())), Color.red);
            
        } else {
            super.mousePressed(e);
            if (sliderDragging) {
                final double val = super.getValue();
                final Color fillColor = (val > System.currentTimeMillis()) ? COLOR_FUTURE : COLOR_PAST;             
                timePopup.show(valueToPosition(val), valueAsString(val), fillColor);
                dragAnim.init(e.getPoint(), e.getWhen());
            } 
        }
    }
    
    /**
     * Deplacement de l'un des sliders
     * @param e 
     */
    @Override
    public void mouseDragged(MouseEvent e) {
        final int x = Math.min(Math.max(0,e.getX()), getWidth());
        if (minLimitDragging) {
            double val = super.positionToValue(x);        
            final DateTime date = SortedListTools.findClosest(availableDates, new DateTime((long)val));
            moveMinLimit(date, true);
            
        } else if (maxLimitDragging) {
            double val = super.positionToValue(x);
            final DateTime date = SortedListTools.findClosest(availableDates, new DateTime((long)val));
            moveMaxLimit(date, true);
        }
        else {
            super.mouseDragged(e);          
                            
            if (sliderDragging) {
                double val = super.getValue();
                timePopup.move(valueToPosition(val), valueAsString(val), dateAsColor((long)val));
                dragAnim.updateVelocity(e.getPoint(), e.getWhen());
            }
        }
    }
    
    /**
     * Sur mouse release, on termine le mouvement de drag si il y a lieu
     * @param e 
     */
    @Override
    public void mouseReleased(MouseEvent e) {
        timePopup.dispose();
        
        if (minLimitDragging) {
            onMinMaxChanged(dateMinAsValue, dateMaxAsValue, false);
            minLimitDragging = false;
        } else if (maxLimitDragging) {
            onMinMaxChanged(dateMinAsValue, dateMaxAsValue, false);
            maxLimitDragging = false;
        } else {
            dragAnim.start(e.getWhen());
        }
    }
    
    /**
     * Sur mouse clicked, on selectionne l'un des slider pour ensuite pouvoir le deplacer
     *   - Soit avec les touches clavier
     *   - Soit via la molette de la souris
     * @param e 
     */
    @Override
    public void mouseClicked(MouseEvent e) {
        this.requestFocus();
        
        // Aucune selection possible en mode veille statique
        if (RefreshMode.STATIC.equals(refreshMode)) {
            return;
        }

        minLimitSelected = maxLimitSelected = sliderSelected = false;
        timePopup.dispose(); // masquage de la derniere popup
        
        final Integer 
                posMinLimit = (dateMinAsValue != null) ? valueToPosition(dateMinAsValue) : null,
                posMaxLimit = (dateMaxAsValue != null) ? valueToPosition(dateMaxAsValue) : null;

        if (!isMinLimitFixed && posMinLimit != null && e.getX() < posMinLimit+5 && e.getX() > posMinLimit - 10) {
            minLimitSelected = true;
            timePopup.show(e.getX(), valueOfAvailableAsString(super.positionToValue(e.getX())), Color.red);
            
        } else if (!isMaxLimitFixed && posMaxLimit != null && e.getX() > posMaxLimit-5 && e.getX() < posMaxLimit + 10) {
            maxLimitSelected = true;
            timePopup.show(e.getX(), valueOfAvailableAsString(super.positionToValue(e.getX())), Color.red);
            
        } else {
            final Integer posVal = /*(dateMinAsValue != null) ?*/ valueToPosition(value);// : null;
            if (e.getX() > posVal-7 && e.getX() < posVal+7) { // Si on clique sur le curseur
                sliderSelected = true;
                final Color fillColor = this.dateAsColor((long)value);        
                timePopup.show(posVal, valueOfAvailableAsString(value), fillColor);
            }
            else {
                super.mouseClicked(e);  
            }
        }
    }
    
    /**
     * Sur Mouse Wheel deplacement du slider selectionné
     * @param e 
     */
    @Override
    public void mouseWheelMoved(MouseWheelEvent e) {
        
        if (minLimitSelected) {
            if (e.getWheelRotation() < 0) {
                System.out.println("first available dates="+availableDates.get(0));
                System.out.println("previous dates = "+ new DateTime((long)(double)dateMinAsValue, DateTimeZone.UTC) + "=>" + SortedListTools.findPrevious(availableDates, new DateTime((long)(double)dateMinAsValue, DateTimeZone.UTC)));
                moveMinLimit(SortedListTools.findPrevious(availableDates, new DateTime((long)(double)dateMinAsValue)), false);       
            } else {
                moveMinLimit(SortedListTools.findNext(availableDates, new DateTime((long)(double)dateMinAsValue)), false);   
            }
            
        } else if (maxLimitSelected) {
            if (e.getWheelRotation() < 0) {
                moveMaxLimit(SortedListTools.findPrevious(availableDates, new DateTime((long)(double)dateMaxAsValue)), false);       
            } else {
                moveMaxLimit(SortedListTools.findNext(availableDates, new DateTime((long)(double)dateMaxAsValue)), false);   
            }     
            
        } else if (sliderSelected) {
            if (e.getWheelRotation() < 0) {
                moveSlider(SortedListTools.findPrevious(validTimesWithoutFiltredFrames, new DateTime((long)value)), false);   
            } else {
                moveSlider(SortedListTools.findNext(validTimesWithoutFiltredFrames, new DateTime((long)value)), false);   
            }            
        }
    }
    
    /**
     * Sur Mouse Move simple mise en lumiere des elements selectionnables (sliders)
     * @param e 
     */
    @Override
    public void mouseMoved(MouseEvent e) {
        final boolean 
                newMinLimitOver,
                newMaxLimitOver;
 
        if (dateMinAsValue != null) {
            final double posMinLimit = valueToPosition(dateMinAsValue);
            newMinLimitOver = (!isMinLimitFixed && e.getX() < posMinLimit+5 && e.getX() > posMinLimit - 10);
        } else {
            newMinLimitOver = false;
        }

        if (dateMaxAsValue != null) {
            final double posMaxLimit = valueToPosition(dateMaxAsValue);
            newMaxLimitOver = (!isMaxLimitFixed && e.getX() > posMaxLimit-5 && e.getX() < posMaxLimit + 10);
        } else {
            newMaxLimitOver = false;
        }
        
        if (newMinLimitOver != minLimitOver || newMaxLimitOver != maxLimitOver) {
            maxLimitOver = newMaxLimitOver;
            minLimitOver = newMinLimitOver;
            repaint();
        }
        
        super.mouseMoved(e);
    }

    /**
     * Sur Mouse Exit, on eneve la surbrillance donnée par le mouse move
     * @param e 
     */
    @Override
    public void mouseExited(MouseEvent e) {
        minLimitOver = maxLimitOver = false;
        repaint();
    }
    
//==============================================================================    
//          +----------------------------------------+
//          |        GESTION DU CLAVIER              |        
//          +----------------------------------------+

    @Override
    public void keyPressed(KeyEvent ke) {
        final int key = ke.getKeyCode();
        if (key == KeyEvent.VK_LEFT) {
            if (minLimitSelected) {
                moveMinLimit(SortedListTools.findPrevious(availableDates, new DateTime((long)(double)dateMinAsValue)), false);       
            } else if (maxLimitSelected) {
                moveMaxLimit(SortedListTools.findPrevious(availableDates, new DateTime((long)(double)dateMaxAsValue)), false);       
            } else if (sliderSelected) {
                moveSlider(SortedListTools.findPrevious(validTimesWithoutFiltredFrames, new DateTime((long)value)), true);   
            }
        } else if (key == KeyEvent.VK_RIGHT) {
            if (minLimitSelected) {
                moveMinLimit(SortedListTools.findNext(availableDates, new DateTime((long)(double)dateMinAsValue)), false); 
            } else if (maxLimitSelected) {
                moveMaxLimit(SortedListTools.findNext(availableDates, new DateTime((long)(double)dateMaxAsValue)), false);       
            } else if (sliderSelected) {
                moveSlider(SortedListTools.findNext(validTimesWithoutFiltredFrames, new DateTime((long)value)), true);   
            }
        } else if (key == KeyEvent.VK_ESCAPE) {
            minLimitSelected = false;
            maxLimitSelected = false;
            sliderSelected = false;
            timePopup.dispose();
        }
    }
    
    /**
     * Gestion de la fin du mouvement 
     * Rq : ne fonctionne pas correctement sous linux repete sans arret
     * => Dans notre cas, ca doublerais toute les requetes au lieu de les allegers
     * @param e 
     */
    @Override public void keyReleased(KeyEvent e) {
//        final int key = e.getKeyCode();
//        
//        if ((key == KeyEvent.VK_LEFT) || (key == KeyEvent.VK_RIGHT)) {
//            if (minLimitSelected) {
//                moveMinLimit(new DateTime((long)(double)dateMinAsValue), false);       
//            } else if (maxLimitSelected) {
//                moveMaxLimit(new DateTime((long)(double)dateMaxAsValue), false);       
//            } else if (sliderSelected) {
//                moveSlider(new DateTime((long)value), false);   
//            }
//        } 
    }

    @Override public void keyTyped(KeyEvent e) {}
           
        
//==============================================================================
//         +----------------------------------------+
//         |         GESTION DU FOCUS               |        
//         +----------------------------------------+

    @Override
    public void focusGained(FocusEvent e) {
    }

    @Override
    public void focusLost(FocusEvent e) {
        minLimitSelected = maxLimitSelected = sliderSelected = false;
        timePopup.dispose();
    }

//==============================================================================
//         +----------------------------------------+
//         |      GESTION DU LANCE DE CURSEUR       |        
//         +----------------------------------------+

    /**
     * Animation permettant de jeter le curseur de position
     */
    private DragAnim dragAnim = new DragAnim(new DragAnim.DragAnimInterface() {
            @Override
            public void onMove(Point2D newPos, double dx, double dy) {
                if (sliderDragging) {
                    setValue(positionToValue((int)newPos.getX()));
                    onDragging();
                    repaint();
                }
            }
            @Override
            public void onEnd() {
                if (sliderDragging) {
                    final Object dragValue = getValueForEvent();
                    if (draggedValueMem != null && !draggedValueMem.equals(dragValue)) {
                        firePropertyChange(VALUE_CHANGED, draggedValueMem, dragValue);
                    }
                    draggedValueMem = draggingValueMem = null;
                    sliderDragging = false;
                }
            }
        });

    @Override
    public Component getComponent() {
        return this;
    }
   
}
