#ifndef CGAL_ARR_CACHED_IO_H
#define CGAL_ARR_CACHED_IO_H

#include <iostream>
#include <tuple>
#include <vector>
#include <type_traits>

namespace CGAL
{

namespace ArrangementIO
{
struct CacheRef
{
  long type;
  long index;

  void set_invalid()
  {
    this->type = -1;
    this->index = -1;
  }

  bool is_invalid() { return this->type == -1 || this->index == -1; }
};

static inline std::ostream& operator<<(std::ostream& os, const CacheRef& cr)
{
  return os << "@(" << cr.type << ' ' << cr.index << ')';
}

static inline std::istream& operator>>(std::istream& is, CacheRef& cr)
{
  while (isspace(is.peek())) is.get();
  if (is.peek() != '@')
  {
    cr.set_invalid();
    return is;
  }

  is.get(); // get '@'
  is.get(); // get '('
  is >> cr.type;
  is >> cr.index;
  is.get(); // get ')'

  return is;
};

template <typename T, typename... Ts>
struct TypeInList : public std::false_type
{
};

template <typename T, typename S, typename... Ts>
struct TypeInList<T, S, Ts...> :
    public std::conditional_t<
      std::is_same<T, S>::value, std::true_type, TypeInList<T, Ts...>>
{
};

template <typename T, typename... Ts>
struct IndexInList
{
  static constexpr int value = 1;
};

template <typename T, typename S, typename... Ts>
struct IndexInList<T, S, Ts...>
{
  static constexpr int value =
    std::is_same<T, S>::value ? 0 : 1 + IndexInList<T, Ts...>::value;
};

template <typename T>
struct TypeIO
{
  template <typename... Ts>
  using void_t = void;

  template <typename S, typename = void>
  struct has_inner_stream : std::false_type
  {
  };

  template <typename S>
  struct has_inner_stream<S, void_t<decltype(&S::inner_stream)>> :
      std::true_type
  {
  };

  template <typename Ostream>
  std::enable_if_t<has_inner_stream<Ostream>::value>
  write(Ostream& os, const T& obj)
  {
    os.inner_stream() << obj;
  }

  template <typename Ostream>
  std::enable_if_t<!has_inner_stream<Ostream>::value>
  write(Ostream& os, const T& obj)
  {
    os << obj;
  }

  template <typename Istream>
  std::enable_if_t<has_inner_stream<Istream>::value> read(Istream& is, T& obj)
  {
    is.inner_stream() >> obj;
  }

  template <typename Istream>
  std::enable_if_t<!has_inner_stream<Istream>::value> read(Istream& is, T& obj)
  {
    is >> obj;
  }
};

template <typename... CachedTypes>
class Cache
{
public:
  template <typename T>
  static constexpr bool is_type_cached()
  {
    return TypeInList<std::decay_t<T>, CachedTypes...>::value;
  }

  template <typename T>
  std::enable_if_t<is_type_cached<T>()> store(const T& obj)
  {
    constexpr int indexInList = IndexInList<T, CachedTypes...>::value;
    auto& container = std::get<indexInList>(storage);
    container.push_back(obj);
  }

  template <typename T>
  std::enable_if_t<!is_type_cached<T>()> store(const T& obj)
  {
  }

  template <typename T>
  std::enable_if_t<is_type_cached<T>(), boost::optional<CacheRef>>
  lookup(const T& obj)
  {
    constexpr int indexInList = IndexInList<T, CachedTypes...>::value;
    auto& container = std::get<indexInList>(storage);
    auto it = std::find(container.begin(), container.end(), obj);
    if (it != container.end())
      return CacheRef{indexInList, std::distance(container.begin(), it)};
    else
      return {};
  }

  template <typename T>
  std::enable_if_t<!is_type_cached<T>(), boost::optional<CacheRef>>
  lookup(const T& obj)
  {
    return {};
  }

  template <typename T>
  std::enable_if_t<is_type_cached<T>(), T> get(CacheRef cache_ref)
  {
    constexpr int indexInList =
      IndexInList<T, CachedTypes...>::value;
    auto& container = std::get<indexInList>(storage);
    return container[cache_ref.index];
  }

private:
  std::tuple<std::vector<CachedTypes>...> storage;
};

template <typename GeomTraits_>
struct TraitsIOCache
{
  using type = Cache<>;
};

template <typename OstreamCache_ = Cache<>>
class CachedOStreamWrapper
{
public:
  using OstreamCache = OstreamCache_;

  CachedOStreamWrapper(std::ostream& os_) : os{&os_} { }
  CachedOStreamWrapper() : os{nullptr} { }

  template <typename T>
  CachedOStreamWrapper& operator<<(const T& obj)
  {
    auto cache_ref = cache.lookup(obj);
    if (cache_ref) { *os << *cache_ref; }
    else
    {
      TypeIO<T>{}.write(*this, obj);
      cache.store(obj);
    }
    return *this;
  }

  decltype(auto) operator<<(std::ostream& (&arg)(std::ostream&))
  {
    *os << arg;
    return *this;
  }

  std::ostream& inner_stream() { return *this->os; }

  operator std::ostream &() const { return *os; }
  operator std::ios &() const { return *os; }
  operator bool() const { return os != nullptr; }
  std::ostream& operator*() const { return *os; }
  std::ostream* operator->() { return os; }

private:
  std::ostream* os;
  OstreamCache cache;
};

template <typename IstreamCache_ = Cache<>>
class CachedIStreamWrapper
{
public:
  using IstreamCache = IstreamCache_;

  CachedIStreamWrapper(std::istream& is_) : is{&is_} { }
  CachedIStreamWrapper() : is{nullptr} { }

  template <typename T>
  std::enable_if_t<
    IstreamCache::template is_type_cached<T>(), CachedIStreamWrapper>&
  operator>>(T& obj)
  {
    CacheRef cache_ref;
    *is >> cache_ref;
    if (cache_ref.is_invalid())
    {
      TypeIO<T>{}.read(*this, obj);
      cache.store(obj);
    }
    else
    {
      obj = this->cache.template get<T>(cache_ref);
    }
    return *this;
  }

  template <typename T>
  std::enable_if_t<
    !IstreamCache::template is_type_cached<T>(), CachedIStreamWrapper>&
  operator>>(T& obj)
  {
    TypeIO<T>{}.read(*this, obj);
    return *this;
  }

  std::istream& inner_stream() { return *this->is; }

  operator std::istream &() const { return *is; }
  operator std::ios &() const { return *is; }
  operator bool() const { return is != nullptr; }
  std::istream& operator*() { return *is; }
  std::istream* operator->() { return is; }

private:
  std::istream* is;
  IstreamCache cache;
};
} // namespace ArrangementIO

} // namespace CGAL

#endif
