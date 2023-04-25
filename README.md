# NEOs_classification

This project consists of a binary classification problem of celestial bodies called NEOs. These are objects belonging to the Solar System whose orbit can intersect the Earth’s orbit and therefore can represent a danger for our planet. NEOs are therefore classified into “Hazardous” and “Not Hazardous”, based on parameters that take into account their potential approach to the Earth, represented by the attributes of the proposed dataset. The aim of this analysis is therefore to define whether these celestial bodies can be dangerous or not for the Earth and to understand which characteristics are useful for this goal.

The proposed dataset consists of 4687 objects whose characteristics are described by 40 attributes. 
These are orbital parameters with mostly real values, suitable for describing the orbit of the body and its structure.

![planetoid-asteroid](https://user-images.githubusercontent.com/81876723/185800413-76597166-4dd7-43fd-a48c-eb468af05e53.jpg)

After the pre-processing part and after some considerations regarding the variables and their relationships, some models were introduced to classify celestial bodies in hazardous and not hazardous.
The models were used with multiple threshold values and both on the original data and on the balanced ones, in order to compare them and evaluate the best one based on the overall error rate and the false negative rate.
