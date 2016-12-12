PIC 10C Final Project
Brian Hsu 

Table of Contents
Background
Program Function
C++ Programming Comments
Demonstration 


Background:
One of the assets most commonly traded on the stock market is called a call option, which gives owners the right 
but not the obligation to purchase a stock at a certain time (the expiration time) for a certain price (the 
strike price). The price of a call option is theoretically supposed to reflect the expected value of the payoff 
(hence the well known Black-Scholes formula), but as with many things, there is often a gap between the theory 
and reality; the purpose of this program is to first quantify the difference, and based on that, roughly 
determine how long of a "lookback period" is most ideal to quantify an option price.


Program Function:
The gist of the program is that it uses real world data (actual prices) to determine some useful constants and 
then uses those constants inside well established formulas from quantitative finance to determine the theoretical 
option price. The process is briefly described below.
1) The function reads in 3 text files, one for the call option prices at different strike prices, one for put 
options, and one for the stock's historical prices  (in the last year). All of this information is then organized 
into a customized matrix object for ease of access and consistency. Next, we use these real world data sets to 
determine some constants.

2) The risk free interest rate is found using the "Put-Call Parity" formula, which basically states that the put 
and call option prices for one strike price should be a certain distance apart to avoid arbitrage opportunities. 
We simply plug in all known variables and solve for the risk free interest rate. 

3) The historical volatility (AKA variance of relative stock price movement) is computed from the historical 
stock data. This is where the user has some choice: the volatility can be extracted from periods of different 
lengths, specifically: week, month, quarter, and year. In a very rough sense, adjusting this lever allows the 
user to determine what lookback period generates the most "true to market" volatility. 

4) The implied volatility is a part of the Black-Scholes formula (which is the most well known theoretical option 
pricing method). Because no closed form of the implied volatility exists, backing it out of the Black-Scholes 
formula requires the real world option price, the previously calculated risk free rate, and performing a Newton-
Raphson approximation.

5) With these useful constants, we can now plug them into the Black-Scholes formula to solve for the theorical 
option price for a user-specified strike price. The formula is computed using both the historical volatility and 
implied volatility, which thereby produces varied results based on what the user selected the lookback period to 
be.

6) Additionally, the binomial tree method (which is essentially the discrete form of the Black-Scholes equation) 
is used to compute the strike price, solely to demonstrate convergence from the discrete to continuous model as 
the size of the timesteps go to zero. 



C++ Programming Comments:
In writing this program, I have endeavored to utilize a variety of tools I've learned from PIC10C this quarter. 
Below I list a few of these programming ideas and how they are used.

1) Inheritance and Class Structure: A custom matrix object is used to read and store much of the data from the 
text files. However, since stocks and options inherently have different parameters, they have seperate classes 
that inherit from the base custom matrix class and include addional details to account for these differences.

2) Templating: Templating is used in many ways, both to generalize operations (such as printing) and for setting 
policy. One instance of this is in finding the risk free rate, the code allows for different kinds of operations 
(max, min, arithmatic average, geometric average) depending on the function pointer that is passed into the 
policy slot of the function.

3) Lambda Functions and Generic Algorithms: These two are used together for almost every vectorized operation in 
the program such as printing and finding averages.

4) A variety of QT tools were implemented to encourage user input and iteraction, which are showcased in the 
demonstration.

Demonstration:
For demonstration of the program, there are three text files which I copied from Google finances and formatted.
Please download the three files and save them in whichever directory the user desires.
Googl_Call_Options.txt
Googl_Put_Options.txt
Google_Stock_History.txt

2) Run the program in QT. Click the "Select ... File" to open the file explorer and select the appropriate file 
for each of the three required files. 

3) Enter the characteristics into the text box (without the quotes)
Stock name: "Google"
Current date: "11.20.2016"
Current price: "760.54"

4) Assemble the matrices by clicking configure

5) Select the strike price and lookback period in the dropdown menus that will populate after appropriately 
loading in all the information

6) Click compute and the results of the computations will be displayed in the large text box on the right!
