# Primary numbers
### University task
##### C++, OpenMP
#### Running
Compile and run using any cpp compiler.
#### Description
The task was to make a programm that calculates Primary Numbers as fast as we can using Sieve of Eratosthenes.

#### Side Notes
Now, the program has one major flaw that I know about (there is no reason for me to fix it since its already finished task)  
Instead of loading chunks of data to memory and working with that data as much as it can be done, it loops throug entire calculation range,  
effectively reloading the same data over and over, creating a huge memory bottleneck, and as such there is no improvement with usage of parallel processing.  
This can be fixed relatively eaisly by processing data in chunks, but I have no reason (and time really) to do this unless its needed for some reason.
