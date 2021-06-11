#ifndef LAB_MC658_PICKUPDELIVERYDECODER_HPP
#define LAB_MC658_PICKUPDELIVERYDECODER_HPP

#include <list>
#include <vector>
#include <algorithm>

class PickupDeliveryDecoder {
public:
  PickupDeliveryDecoder();
  ~PickupDeliveryDecoder();

  double decode(const std::vector<double> &chromosome) const;
};

#endif // LAB_MC658_PICKUPDELIVERYDECODER_HPP
