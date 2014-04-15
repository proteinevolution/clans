package clans.gui;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.Timer;
import javax.swing.UIManager;

public class GuiMessageOverlay extends JComponent {
	private static final long serialVersionUID = 818220626180736376L;

	/**
	 * States of the message overlay, e.g. OFF, LOADING, SAVING_COMPLETED.
	 */
	enum State {
		OFF,
		CUSTOM,
		
		LOADING,
		LOADING_COMPLETED,
		LOADING_CANCELED,
		LOADING_FAILED,
		
		SAVING,
		SAVING_COMPLETED,
		SAVING_CANCELED,
		SAVING_FAILED;
	}

	/**
	 * Messages for different {@code State}s.
	 */
	enum Message {
		OFF("MESSAGE OVERLAY OFF"),
		
		CUSTOM("CUSTOM MESSAGE SET"),
		
		LOADING("loading file, please wait"),
		LOADING_COMPLETED("loaded file successfully"),
		LOADING_CANCELED("canceled loading file"),
		LOADING_FAILED("LOADING FILE FAILED!"),

		SAVING("saving file, please wait"),
		SAVING_COMPLETED("saved file successfully"),
		SAVING_CANCELED("canceled saving file"),
		SAVING_FAILED("SAVING FILE FAILED!");

		private String value;

		private Message(String value) {
			this.value = value;
		}
		
		protected String get() {
			return value;
		}
	}
	
	enum MessageDetails {
		LOADING("press escape to cancel"),
		SAVING("press escape to cancel");

		private String value;

		private MessageDetails(String value) {
			this.value = value;
		}
		
		protected String get() {
			return value;
		}
	}

	/**
	 * Durations in milliseconds for different message types, e.g. INFO or ERROR.
	 */
	enum Duration {
		
		FADING(2000),
		INFO(2500),
		WARNING(3500),
		ERROR(4500),
		INFINITE(Integer.MAX_VALUE);
		
		private int value;
		
		private Duration(int value) {
			this.value = value;
		}
		
		protected int get() {
			return value;
		}
	}
	

	enum TextColor {
		WORKING(255, 128, 0), // an orange
		SUCCESS(128, 255, 0), // a green
		ERROR(160, 0, 0); // a dark red

		private Color value;

		private TextColor(int red, int green, int blue) {
			this.value = new Color(red, green, blue);
		}
		
		protected Color get() {
			return value;
		}
	}

	private Component top_spacer; // moves the labels down
	private JLabel mainLabel; // shows the main message
	private JLabel detailsLabel; // shows messageDetails
	private int labelInsetTop; // used for custom spacing of the labels
	private int labelInsetBottom; // used for custom spacing of the labels
	private int labelInsetLeftRight; // used for custom spacing of the labels
	private Color colorLabelBackground; // the labels background is not translucent for better visibility

	private Timer timer; // draw only every X milliseconds instead of using up one CPU for drawing continuously.


	private long startTime; // time of last state change
	private State state; // the current state
	private String currentMessage; // the details of the message, e.g. what went wrong during loading
	private int currentDuration; // the duration for the current mode
	private boolean currentIsFading; // whether the current message fades or abruptly vanishes
	private int currentFading; // the fading duration
	private boolean currentHasProgressDots; // whether the current message gets appended progress-indicating dots

	private String currentDots; // String representation of the currently required progress-indicating dots (can be "")
	private int dotTicks; // the number of ticks seen since the last update of currentDots
	final private int millisecondsPerDot = 350; // milliseconds that must pass for the next different progress display
	
	private float alpha; // the transparency of the message overlay
	
	public GuiMessageOverlay() {

		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		colorLabelBackground = new Color(128, 128, 128);
		
		labelInsetTop = 4;
		labelInsetBottom = 2;
		labelInsetLeftRight = 4;
				
		mainLabel = new JLabel();
		mainLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(20f));
		mainLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		mainLabel.setOpaque(true);
		mainLabel.setBackground(colorLabelBackground);
		mainLabel.setBorder(BorderFactory.createEmptyBorder(labelInsetTop, labelInsetLeftRight, labelInsetBottom,
				labelInsetLeftRight));
		
		detailsLabel = new JLabel();
		detailsLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(16f));
		detailsLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		detailsLabel.setBorder(BorderFactory.createEmptyBorder(labelInsetTop, labelInsetLeftRight, labelInsetBottom,
				labelInsetLeftRight));
		detailsLabel.setOpaque(true);
		detailsLabel.setBackground(colorLabelBackground);

		// the initial values for top_spacer ARE NOT USED as we overwrite them in updateGlassSize() below
		top_spacer = Box.createRigidArea(new Dimension(1, 1));
		add(top_spacer);
		add(mainLabel);
		add(detailsLabel);

		currentFading = Duration.FADING.get(); // this is kept stable for now for consistent fading looks

		int timer_delay = 100; // in milliseconds
		setupTimer(timer_delay);
	
		setOff();
	}

	protected void setState(State state) {
		this.state = state;
		startTime = System.currentTimeMillis();
	}

	/**
	 * Sets up Timer to only redraw every X milliseconds instead of continuously, which would use up 1 CPU entirely just
	 * for the overlay.
	 */
	private void setupTimer(int delay) {
		timer = new Timer(delay, new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				invoke();
			}
		});
		timer.setInitialDelay(0);
	}

	/**
	 * Checks if and which type of message is waiting to be displayed and initiates its display. Does nothing it no
	 * messages is about to be displayed. Disables the message overlay if message duration is over.
	 */
	private void invoke() {
		
		if (state == State.OFF) {
			return;
		}
		
		// time since start is also relevant for fading, hence we forward it for faded message (see below)
		long time_since_start = System.currentTimeMillis() - startTime;
		if (time_since_start > currentDuration) {
			setOff();
			return;
		}

		resetTransparency();

		if (currentIsFading) {
			showFadingMessage(time_since_start);
		} else {
			showMessage();
		}
	}
	
	/**
	 * Sets the main message.
	 * @param message The main message text.
	 */
	private void setMessage(String message){
		mainLabel.setText(message);
		currentMessage = message;
	}
	
	/**
	 * Resets message to none.
	 */
	private void resetMessage() {
		setMessage("");
	}
	
	/**
	 * Sets message details.
	 * 
	 * @param details
	 *            The details message.
	 */
	protected void setMessageDetails(String details) {

		if (details.length() > 70) { // shorten long messages
			details = details.substring(0, 65) + "[...]";
		}
		
		detailsLabel.setText(details);
		detailsLabel.setVisible(true);
	}
	
	/**
	 * Resets detailed message to none.
	 */
	private void resetMessageDetails() {
		detailsLabel.setVisible(false);
	}

	/**
	 * Sets the text color for main and details message.
	 * @param color The new color.
	 */
	private void setColor(Color color) {
		mainLabel.setForeground(color);
		detailsLabel.setForeground(color);
	}
	
	/**
	 * Resets text color to default (currently: a green)
	 */
	private void resetColor() {
		setColor(TextColor.SUCCESS.get());
	}
	
	/**
	 * Sets the duration of the current message.
	 * 
	 * @param duration
	 *            The message duration.
	 */
	private void setDuration(Duration duration) {
		currentDuration = duration.get();
	}
	
	/**
	 * Resets the duration to infinite.
	 * 
	 */
	private void resetDuration() {
		setDuration(Duration.INFINITE);
	}
	
	/**
	 * Enables fading for the message.
	 */
	private void enableFading() {
		currentIsFading = true;
	}
	
	/**
	 * Disables fading for the message.
	 */
	private void disableFading() {
		currentIsFading = false;
	}
	
	/**
	 * Enables progress-indicating dots for this message.   
	 */
	private void enableProgressDots() {
		currentDots = "";
		currentHasProgressDots= true;
	}
	
	/**
	 * Disables progress-indicating dots for this message.   
	 */
	private void disableProgressDots() {
		currentDots = "";
		currentHasProgressDots = false;
	}

	/**
	 * Resets the transparency of the overlay to "not transparent at all".
	 */
	private void resetTransparency() {
		alpha = 1;
		repaint();
	}
	
	/**
	 * Sets the correct number of dots for recurring calls to indicate work in progress.
	 * 
	 * @return true if the number of dots changes and hence the main label must be updated, false else.
	 */
	private boolean updateDots() {
		dotTicks += 1;

		if (dotTicks * timer.getDelay() > millisecondsPerDot) {
			if (currentDots.length() < 3) {
				currentDots += ".";
			} else {
				currentDots = "";
			}

			dotTicks = 0;
			return true;
		
		} else{
			return false;
		}
	}
	
	/**
	 * Switches the overlay off for the time being. Used once in constructor and whenever a message does not need to be
	 * displayed any more.
	 */
	private void setOff() {
		setVisible(false);
		
		if (timer != null) {
			timer.stop();
		}

		resetMessage();
		resetMessageDetails();
		resetColor();
		resetDuration();
		
		disableFading();
		disableProgressDots();
		resetTransparency();
		
		setState(State.OFF);
	}

	/**
	 * Changes the state to "loading has completed". Can contain a detailed message.
	 */
	private void setLoadingCompleted(String details) {
		setupStandardMessage(State.LOADING_COMPLETED, details);
	}

	/**
	 * Changes the state to "loading was canceled". Can contain a detailed message.
	 */
	private void setLoadingCanceled(String details) {
		setupStandardMessage(State.LOADING_CANCELED, details);
	}

	/**
	 * Changes the state to "loading has failed". Can contain a detailed message.
	 */
	private void setLoadingFailed(String details) {
		setupStandardMessage(State.LOADING_FAILED, details);
	}

	/**
	 * Changes the state to "saving has completed". Can contain a detailed message.
	 */
	private void setSavingCompleted(String details) {
		setupStandardMessage(State.SAVING_COMPLETED, details);
	}

	/**
	 * Changes the state to "saving was canceled". Can contain a detailed message.
	 */
	private void setSavingCanceled(String details) {
		setupStandardMessage(State.SAVING_CANCELED, details);
	}

	/**
	 * Changes the state to "saving has failed". Can contain a detailed message.
	 */
	private void setSavingFailed(String details) {
		setupStandardMessage(State.SAVING_FAILED, details);
	}

	/**
	 * Displays the message. If applicable, the correct number of progress-indicating dots are set here.
	 */
	protected void showMessage() {
	
		if (currentHasProgressDots && updateDots()) {
			mainLabel.setText(currentMessage + currentDots);
		}
		
		setVisible(true);
	}

	/**
	 * Displays the message (and the message details, if any) and fades the color according to the fraction of time of
	 * the complete duration for the message.
	 * 
	 * @param time_since_start
	 *            Time that has passes since this message was first shown.
	 */
	private void showFadingMessage(long time_since_start) {
		// fade by reducing visibility for the whole glass pane; see paintComponent below.
		alpha = 1 - getFadingValue(time_since_start);
		
		showMessage();
		repaint();
	}

	/**
	 * Computes the fading value. After an initial delay of (Message duration - fading duration) seconds in which no
	 * fading (==0) is returned, the fraction of time already elapsed w.r.t. the complete message duration is returned.
	 * The speed of fading from 0 -> 1 values is linear.
	 * 
	 * @param time_since_start
	 * @return Fraction of duration already passed as value in [0, 1]
	 */
	private float getFadingValue(long time_since_start) {
		
		float timespan_without_fading = currentDuration - currentFading;
		
		// make the initial delay period non-transparent
		if (time_since_start < timespan_without_fading) {
			return 0;
		}

		// fade during the last part the duration and make the subtraction negative number proof
		float fade_fraction = (float) Math.max(time_since_start - timespan_without_fading, 0) / currentFading;
		return fade_fraction;
	}

	/**
	 * Adjusts the label placement when the overlay's parent is resized by changing the top spacer height. Used once in
	 * constructor.
	 * 
	 * @param new_height
	 *            The new height of the overlay
	 */
	protected void updateGlassSize(int new_height) {
		// top_spacer is the first entry
		remove(0);

		top_spacer = Box.createRigidArea(new Dimension(0, new_height / 3 * 2));

		add(top_spacer, 0);

		validate();
	}


	/**
	 * Shows message if there is one, otherwise does nothing.
	 */
	protected void activate_overlay() {
		if (state != State.OFF && !timer.isRunning()) {
			timer.start();
			
		}
	}

	/**
	 * Sets up a message with all available details of the overlay.
	 * 
	 * @param state The state of the overlay for this new message.
	 * @param main_message The main message text.
	 * @param message_details The message details text.
	 * @param color The text color.
	 * @param duration The message duration.
	 * @param is_fading true if the message is supposed to fade, false else
	 * @param has_progress_dots true if the message is supposed to have progress-indicating dots.
	 */
	private void setupMessage(State state, String main_message, String message_details, Color color,
			Duration duration, boolean is_fading, boolean has_progress_dots) {

		setState(state);
		
		setMessage(main_message);		
		
		if (message_details == null) {
			resetMessageDetails();
		} else {
			setMessageDetails(message_details);
		}
		
		setColor(color);
		
		setDuration(duration);

		if (is_fading) {
			enableFading();
		} else {
			disableFading();
		}
		
		if (has_progress_dots) {
			enableProgressDots();
		} else {
			disableProgressDots();
		}
		
		activate_overlay();
	}
	
	/**
	 * Sets up a message with progress-indicating dots and no fading.
	 * 
	 * @param state
	 *            The state of the overlay for this new message.
	 * @param main_message
	 *            The main message text.
	 * @param message_details
	 *            The message details text.
	 */
	protected void setupProgressMessage(State state, String main_message, String message_details) {
		setupMessage(state, main_message, message_details, TextColor.WORKING.get(), Duration.INFINITE, false, true);
	}
	
	/**
	 * Sets up a message with fading and no progress-indicating dots.
	 * 
	 * @param state
	 *            The state of the overlay for this new message.
	 * @param main_message
	 *            The main message text.
	 * @param message_details
	 *            The message details text.
	 * @param color
	 *            The text color.
	 * @param duration
	 *            The message duration.
	 */
	private void setupFadingMessage(State state, String main_message, String message_details, Color color,
			Duration duration) {
		setupMessage(state, main_message, message_details, color, duration, true, false);
	}

	/**
	 * Sets up the standard message according to the given state.
	 * 
	 * @param state
	 *            The state for which the standard message should be shown.
	 */
	private void setupStandardMessage(State state) {
		setupStandardMessage(state, null);
	}
	
	/**
	 * Sets up the standard message according to the given state with custom details message.
	 * 
	 * @param state
	 *            The state for which the standard message should be shown.
	 */
	private void setupStandardMessage(State state, String details) {
		String main_message = Message.valueOf(state.toString()).get();

		switch (state) {
		case OFF:
			return;
			
		case CUSTOM:
			// custom message are set up by their own method
			return;
			
		case LOADING:
		case SAVING:
			setupProgressMessage(state, main_message, MessageDetails.valueOf(state.toString()).get());
			return;
			
		case LOADING_CANCELED:
		case SAVING_CANCELED:
			setupFadingMessage(state, main_message, details, TextColor.WORKING.get(), Duration.INFO);
			return;
			
		case LOADING_COMPLETED:
		case SAVING_COMPLETED:
			setupFadingMessage(state, main_message, details, TextColor.SUCCESS.get(), Duration.INFO);
			return;

		case LOADING_FAILED:
		case SAVING_FAILED:
			setupFadingMessage(state, main_message, details, TextColor.ERROR.get(), Duration.ERROR);
			return;

		default:
			break;
		}
	}

	/**
	 * We change the transparency of the whole overlay during repaints to simulate message fading.
	 */
	protected void paintComponent(Graphics g) {
		Graphics2D g2d = (Graphics2D) g;
		g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, alpha));
	}
	
	/**
	 * @return Default color for completed processes.
	 */
	protected Color getColorSuccess() {
		return TextColor.SUCCESS.get();
	}

	/**
	 * @return The Overlays default color.
	 */
	protected Color getColorWorking() {
		return TextColor.WORKING.get();
	}

	/**
	 * @return Default color for erroneous processes.
	 */
	protected Color getColorError() {
		return TextColor.ERROR.get();
	}
	
	/**
	 * @return Default duration for information-type messages.
	 */
	protected Duration getDurationInfo() {
		return Duration.INFO;
	}
	
	/**
	 * @return Default duration for information-type messages.
	 */
	protected Duration getDurationWarning() {
		return Duration.WARNING;
	}
	
	/**
	 * @return Default duration for information-type messages.
	 */
	protected Duration getDurationError() {
		return Duration.ERROR;
	}
	
	/**
	 * Sets up a custom message with all available details of the overlay.
	 * 
	 * @param main_message
	 *            The main message text.
	 * @param message_details
	 *            The message details text.
	 * @param color
	 *            The text color.
	 * @param duration
	 *            The message duration.
	 * @param is_fading
	 *            true if the message is supposed to fade, false else
	 * @param has_progress_dots
	 *            true if the message is supposed to have progress-indicating dots.
	 */
	protected void setCustomMessage(String main_message, String message_details, Color color, Duration duration,
			boolean is_fading, boolean with_progress_dots) {
		setupMessage(State.CUSTOM, main_message, message_details, color, duration, is_fading, with_progress_dots);
	}
	
	/**
	 * Changes the state to "currently loading".
	 */
	protected void setLoading() {
		setupStandardMessage(State.LOADING);
	}

	/**
	 * Changes the state to "currently saving".
	 */
	protected void setSaving() {
		setupStandardMessage(State.SAVING);
	}

	/**
	 * Changes the state to "loading/saving has completed" depending on previous state.
	 */
	protected void setCompleted() {
		setCompleted(null);
	}
	
	/**
	 * Changes the state to "loading/saving has completed" depending on previous state and adds a detailed message.
	 */
	protected void setCompleted(String details) {

		switch (state) {
		case LOADING:
			setLoadingCompleted(details);
			break;

		case SAVING:
			setSavingCompleted(details);
			break;

		default:
			System.err.println("GuiMessageOverlay.setCompleted not implemented while in state \"" + state.toString()
					+ "\".");
			break;
		}
	}

	/**
	 * Changes the state to "loading/saving was canceled" depending on previous state.
	 */
	protected void setCanceled() {
		setCanceled(null);
	}
	
	/**
	 * Changes the state to "loading/saving was canceled" depending on previous state and adds a detailed message.
	 */
	protected void setCanceled(String details) {

		switch (state) {
		case LOADING:
			setLoadingCanceled(details);
			break;

		case SAVING:
			setSavingCanceled(details);
			break;

		default:
			System.err.println("GuiMessageOverlay.setCanceled not implemented while in state \"" + state.toString()
					+ "\".");
			break;
		}
	}

	/**
	 * Changes the state to "loading/saving has failed" depending on previous state.
	 */
	protected void setFailed() {
		setFailed(null);
	}

	/**
	 * Changes the state to "loading/saving has failed" depending on previous state and adds a detailed message.
	 * 
	 * @param details
	 *            A message stating details, e.g. file name or reason for failure.
	 */
	protected void setFailed(String details) {
	
		switch (state) {
		case LOADING:
			setLoadingFailed(details);
			break;

		case SAVING:
			setSavingFailed(details);
			break;

		default:
			System.err.println("GuiMessageOverlay.setFailed not implemented while in state \"" + state.toString()
					+ "\".");
			break;
		}
	}
}

/**
 * NEVER USE THIS CLASS IN A PRODUCTION RELEASE!!!
 * <p>
 * This class wraps {@code GuiMessageOverlay} for debugging and shows various information in stderr.
 */
class GuiMessageOverlayLogged extends GuiMessageOverlay {
	private static final long serialVersionUID = -8396089731118502477L;

	private int message_count;
	private long last_message_time;
	
	public GuiMessageOverlayLogged() {
		super();
		message_count = 0;
		last_message_time = 0;
	}
	
	/**
	 * Shows state changes.
	 */
	protected void setState(State state) {
		super.setState(state);
		System.err.println("state is now: " + Message.valueOf(state.toString()));
		message_count = 0;
	}
	
	/**
	 * Shows statistic on message frequency.
	 */
	protected void showMessage() {
		super.showMessage();
		
		long current_time = System.currentTimeMillis();
		System.err.println("\ttime since last message: " + (current_time - last_message_time));
		last_message_time = current_time;
		
		message_count++;
		System.err.println("\tmessage refresh No. " + message_count);
	}
	
	protected void setupProgressMessage(State state, String main_message, String message_details) {
		System.err.println("setupProgressmessage(" + state.toString() + ", \"" + main_message + "\", \""
				+ message_details + "\")");
		super.setupProgressMessage(state, main_message, message_details);
	}
}
